#include "object_tracker.h"

#include <logger.h>

#include "common/geodetic/geodetic_coordinates.h"
#include "common/tracking/tle.h"
#include "core/plugin.h"
#include "common/utils.h"

namespace satdump
{
    ObjectTracker::ObjectTracker(bool is_gui) : is_gui(is_gui)
    {
        if (general_tle_registry.size() > 0)
            has_tle = true;

        for (auto &tle : general_tle_registry)
            satoptions.push_back(tle.name);

        satellite_observer_station = predict_create_observer("Main", 0, 0, 0);

        // Updates on registry updates
        eventBus->register_handler<TLEsUpdatedEvent>([this](TLEsUpdatedEvent)
                                                     {
                                                            general_mutex.lock();

                                                            if (general_tle_registry.size() > 0)
                                                                has_tle = true;

                                                            satoptions.clear();
                                                            for (auto &tle : general_tle_registry)
                                                                satoptions.push_back(tle.name);
                                                                
                                                            general_mutex.unlock(); });

        // Start threads
        backend_thread = std::thread(&ObjectTracker::backend_run, this);
        rotatorth_thread = std::thread(&ObjectTracker::rotatorth_run, this);
    }

    ObjectTracker::~ObjectTracker()
    {
        backend_should_run = false;
        if (backend_thread.joinable())
            backend_thread.join();

        rotatorth_should_run = false;
        if (rotatorth_thread.joinable())
            rotatorth_thread.join();

        predict_destroy_observer(satellite_observer_station);

        if (satellite_object != nullptr)
            predict_destroy_orbital_elements(satellite_object);


        for (auto & [norad, sat_obj] : norad_to_sat_object)
        {
            predict_destroy_orbital_elements(sat_obj);
        }
    }

    void ObjectTracker::setQTH(double qth_lon, double qth_lat, double qth_alt)
    {
        general_mutex.lock();
        this->qth_lon = qth_lon;
        this->qth_lat = qth_lat;
        this->qth_alt = qth_alt;
        if (satellite_observer_station != nullptr)
            predict_destroy_observer(satellite_observer_station);
        satellite_observer_station = predict_create_observer("Main", qth_lat * DEG_TO_RAD, qth_lon * DEG_TO_RAD, qth_alt);
        backend_needs_update = true;
        general_mutex.unlock();
    }

    void ObjectTracker::setObject(TrackingMode mode, int objid)
    {
        general_mutex.lock();
        tracking_mode = TRACKING_NONE;

        if (mode == TRACKING_HORIZONS)
        {
            if (horizonsoptions.size() == 1)
                horizonsoptions = pullHorizonsList();
            for (int i = 0; i < (int)horizonsoptions.size(); i++)
            {
                if (horizonsoptions[i].first == objid)
                {
                    tracking_mode = TRACKING_HORIZONS;
                    current_horizons_id = i;
                    break;
                }
            }
        }
        else if (mode == TRACKING_SATELLITE)
        {
            for (int i = 0; i < (int)satoptions.size(); i++)
            {
                if (general_tle_registry[i].norad == objid)
                {
                    tracking_mode = TRACKING_SATELLITE;
                    current_satellite_id = i;
                    break;
                }
            }
        }

        backend_needs_update = true;
        general_mutex.unlock();
    }

    void ObjectTracker::setObjects(TrackingMode mode, const std::vector<SatellitePass> &sat_passes)
    {// calculate pass points for each satellite in schedule
        if (mode == TRACKING_SATELLITE)
        {
            upcoming_satellite_passes_mtx.lock();
            upcoming_satellite_passes = sat_passes;
            upcoming_satellite_pass_points.clear(); // clear previous pass points

            // clear previous sat objects
            for (auto & [norad, sat_obj] : norad_to_sat_object)
            {
                predict_destroy_orbital_elements(sat_obj);
            }
            norad_to_sat_object.clear();

            // find all relevant sat objects in registry and update
            for(const SatellitePass& satellite_pass : upcoming_satellite_passes)
            {   // skip already found objects
                if (!norad_to_sat_object.count(satellite_pass.norad))
                {
                    for (auto & tle : general_tle_registry)
                    {
                        if (tle.norad == satellite_pass.norad)
                        {
                            norad_to_sat_object.insert({tle.norad, predict_parse_tle(tle.line1.c_str(), tle.line2.c_str())});
                            break;
                        }
                    }
                }
            }
            // Calculate pass point for each upcoming satellite
            for(const auto& [norad, next_aos_time, next_los_time, max_elevation] : upcoming_satellite_passes)
            {
                const auto satellite_object = norad_to_sat_object[norad];

                std::vector<SatAzEl> curr_sat_upcoming_pass_points;
                if (!predict_is_geosynchronous(satellite_object)) // TODO test on geosynchronous
                {
                    constexpr int time_steps = 50;
                    // Calculate a few points during the pass
                    predict_position satellite_orbit2;
                    predict_observation observation_pos2;

                    double time_step = abs(next_los_time - next_aos_time) / time_steps;
                    for (double ctime = next_aos_time; ctime <= next_los_time; ctime += time_step)
                    {
                        predict_orbit(satellite_object, &satellite_orbit2, predict_to_julian_double(ctime));
                        predict_observe_orbit(satellite_observer_station, &satellite_orbit2, &observation_pos2);
                        curr_sat_upcoming_pass_points.push_back({
                            float(observation_pos2.azimuth * RAD_TO_DEG), float(observation_pos2.elevation * RAD_TO_DEG)
                        });
                    }
                } else
                {
                    // FIXME push single point for geosynchronous?
                    logger->warn("Skipping geosynchronous satellite");
                }
                upcoming_satellite_pass_points.push_back(curr_sat_upcoming_pass_points);
            }
            upcoming_satellite_passes_mtx.unlock();
        } else { logger->warn("Only TRACKING_SATELLITE mode is supported!"); }
    }

    void ObjectTracker::setRotator(std::shared_ptr<rotator::RotatorHandler> rot)
    {
        rotator_handler_mtx.lock();
        rotator_handler = rot;
        rotator_handler_mtx.unlock();
    }

    void ObjectTracker::setRotatorEngaged(bool v)
    {
        rotator_handler_mtx.lock();
        rotator_engaged = v;
        rotator_handler_mtx.unlock();
    }

    void ObjectTracker::setRotatorTracking(bool v)
    {
        rotator_handler_mtx.lock();
        rotator_tracking = v;
        rotator_handler_mtx.unlock();
    }

    void ObjectTracker::setRotatorReqPos(float az, float el)
    {
        rotator_handler_mtx.lock();
        rot_current_req_pos.az = az;
        rot_current_req_pos.el = el;
        rotator_handler_mtx.unlock();
    }
}
