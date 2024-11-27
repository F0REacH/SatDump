#include <logger.h>

#include "object_tracker.h"
#include "common/utils.h"

namespace satdump
{
    nlohmann::json ObjectTracker::getStatus()
    {
        nlohmann::json v;

        std::string obj_name = "None";
        if (tracking_mode == TRACKING_HORIZONS)
            obj_name = horizonsoptions[current_horizons_id].second;
        else if (tracking_mode == TRACKING_SATELLITE)
            obj_name = satoptions[current_satellite_id];
        v["name"] = obj_name;
        v["current_position"] = sat_current_pos;

        if (tracking_mode == TRACKING_SATELLITE && satellite_object != nullptr)
        {
            v["sat_azimuth_rate"] = satellite_observation_pos.azimuth_rate * RAD_TO_DEG;
            v["sat_elevation_rate"] = satellite_observation_pos.elevation_rate * RAD_TO_DEG;
            v["sat_current_range"] = satellite_observation_pos.range;
            v["pass_points"] = getPassPoints(satellite_object, next_aos_time, next_los_time);
        }

        v["aos_time"] = next_aos_time;
        v["los_time"] = next_los_time;

        double ctime = getTime();
        v["next_event_in"] = getNextEventIn(ctime, next_aos_time, next_los_time);
        v["next_event_is_aos"] = next_aos_time > ctime;

        v["rotator_engaged"] = rotator_engaged;
        v["rotator_tracking"] = rotator_tracking;
        v["rot_current_pos"] = rot_current_pos;
        v["rot_current_req_pos"] = rot_current_req_pos;

        return v;
    }

    image::Image ObjectTracker::getPolarPlotImg(int plot_size)
    {
        image::Image img(8, plot_size, plot_size, 3);

        // All black bg
        img.fill(0);

        // Draw the "target-like" plot with elevation rings
        float radius = 0.45;
        float radius1 = plot_size * radius * (3.0 / 9.0);
        float radius2 = plot_size * radius * (6.0 / 9.0);
        float radius3 = plot_size * radius * (9.0 / 9.0);

        std::vector<double> color_green = {0, 1, 0};
        std::vector<double> color_red = {1, 0, 0};
        std::vector<double> color_orange = {1, 165.0 / 255.0, 0};
        std::vector<double> color_cyan = {0, 237.0 / 255.0, 1};

        img.draw_circle(plot_size / 2, plot_size / 2,
                        radius1, color_green, false);
        img.draw_circle(plot_size / 2, plot_size / 2,
                        radius2, color_green, false);
        img.draw_circle(plot_size / 2, plot_size / 2,
                        radius3, color_green, false);

        img.draw_line(plot_size / 2, 0,
                      plot_size / 2, plot_size - 1,
                      color_green);
        img.draw_line(0, plot_size / 2,
                      plot_size - 1, plot_size / 2,
                      color_green);

        // Draw the satellite's trace
        if (upcoming_pass_points.size() > 1)
        {
            upcoming_passes_mtx.lock();
            for (int i = 0; i < (int)upcoming_pass_points.size() - 1; i++)
            {
                auto &p1 = upcoming_pass_points[i];
                auto &p2 = upcoming_pass_points[i + 1];

                float point_x1, point_x2, point_y1, point_y2;
                point_x1 = point_x2 = plot_size / 2;
                point_y1 = point_y2 = plot_size / 2;

                point_x1 += az_el_to_plot_x(plot_size, radius, p1.az, p1.el);
                point_y1 -= az_el_to_plot_y(plot_size, radius, p1.az, p1.el);

                point_x2 += az_el_to_plot_x(plot_size, radius, p2.az, p2.el);
                point_y2 -= az_el_to_plot_y(plot_size, radius, p2.az, p2.el);

                img.draw_line(point_x1, point_y1,
                              point_x2, point_y2,
                              color_orange);
            }
            upcoming_passes_mtx.unlock();
        }

        // Draw the current satellite position
        if (sat_current_pos.el > 0)
        {
            float point_x = plot_size / 2;
            float point_y = plot_size / 2;

            point_x += az_el_to_plot_x(plot_size, radius, sat_current_pos.az, sat_current_pos.el);
            point_y -= az_el_to_plot_y(plot_size, radius, sat_current_pos.az, sat_current_pos.el);

            img.draw_circle(point_x, point_y, 5, color_red, true);
        }

        if (rotator_handler && rotator_handler->is_connected())
        {
            {
                float point_x = plot_size / 2;
                float point_y = plot_size / 2;

                point_x += az_el_to_plot_x(plot_size, radius, rot_current_pos.az, rot_current_pos.el);
                point_y -= az_el_to_plot_y(plot_size, radius, rot_current_pos.az, rot_current_pos.el);

                img.draw_circle(point_x, point_y, 9, color_cyan, false);
            }

            if (rotator_engaged)
            {
                float point_x = plot_size / 2;
                float point_y = plot_size / 2;

                point_x += az_el_to_plot_x(plot_size, radius, rot_current_req_pos.az, rot_current_req_pos.el);
                point_y -= az_el_to_plot_y(plot_size, radius, rot_current_req_pos.az, rot_current_req_pos.el);

                img.draw_line(point_x - 5, point_y, point_x - 12, point_y, color_cyan);
                img.draw_line(point_x + 5, point_y, point_x + 12, point_y, color_cyan);
                img.draw_line(point_x, point_y - 5, point_x, point_y - 12, color_cyan);
                img.draw_line(point_x, point_y + 5, point_x, point_y + 12, color_cyan);
            }
        }

        return img;
    }

    ObjectTracker::SatAzEl ObjectTracker::getPredictionPoint(const predict_orbital_elements_t *satellite_object, const double prediction_time) const
    {
        predict_position satellite_orbit2;
        predict_observation observation_pos2;

        predict_orbit(satellite_object, &satellite_orbit2, predict_to_julian_double(prediction_time));
        predict_observe_orbit(satellite_observer_station, &satellite_orbit2, &observation_pos2);

        return {
            float(observation_pos2.azimuth * RAD_TO_DEG),
            float(observation_pos2.elevation * RAD_TO_DEG)
        };
    }

    std::vector<ObjectTracker::SatAzEl> ObjectTracker::getPassPoints(const predict_orbital_elements_t *satellite_object, const double aos_time, const double los_time, const int time_steps) const
    {
        std::vector<SatAzEl> pass_points;
        if (!predict_is_geosynchronous(satellite_object))
        {
            double time_step = abs(los_time - aos_time) / time_steps;
            for (double ctime = aos_time; ctime <= los_time; ctime += time_step)
            {
                pass_points.push_back(getPredictionPoint(satellite_object, ctime));
            }
        }
        else
        {
            logger->warn("Skipping geosynchronous satellite");
        }
        return pass_points;
    }

    double ObjectTracker::getNextEventIn(double ctime, double aos_time, double los_time)
    {
        double timeOffset = 0;
        if (aos_time > ctime)
            timeOffset = aos_time - ctime;
        else
            timeOffset = los_time - ctime;
        return timeOffset;
    }

    nlohmann::json ObjectTracker::getTrackedSatelliteObjects()
    {
        tracked_satellite_objects_mtx.lock();
        // update satellite positions if needed

        for(auto& tr_obj : tracked_satellite_objects)
        {
            if (norad_to_sat_object.count(tr_obj.norad))
            {
                const double current_time = getTime();
                tr_obj.next_event_in = getNextEventIn(current_time, tr_obj.aos_time, tr_obj.los_time);
                tr_obj.next_event_is_aos = tr_obj.aos_time > current_time;

                // calculate only within time window
                if (current_time > tr_obj.aos_time && current_time < tr_obj.los_time)
                {
                    tr_obj.current_position = getPredictionPoint(norad_to_sat_object[tr_obj.norad], current_time);
                } else
                {
                    tr_obj.current_position = SatAzEl();
                }
            } else
            {
                logger->warn("Could not find norad in registry! norad:"+std::to_string(tr_obj.norad));
            }
        }
        nlohmann::json v = tracked_satellite_objects;
        tracked_satellite_objects_mtx.unlock();
        return v;
    }
}
