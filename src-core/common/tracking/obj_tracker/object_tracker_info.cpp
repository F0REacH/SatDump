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
        v["object_name"] = obj_name;
        v["sat_current_pos"] = sat_current_pos;

        if (tracking_mode == TRACKING_SATELLITE && satellite_object != nullptr)
        {
            v["sat_azimuth_rate"] = satellite_observation_pos.azimuth_rate * RAD_TO_DEG;
            v["sat_elevation_rate"] = satellite_observation_pos.elevation_rate * RAD_TO_DEG;
            v["sat_current_range"] = satellite_observation_pos.range;
        }

        v["next_aos_time"] = next_aos_time;
        v["next_los_time"] = next_los_time;

        double timeOffset = 0, ctime = getTime();
        if (next_aos_time > ctime)
            timeOffset = next_aos_time - ctime;
        else
            timeOffset = next_los_time - ctime;

        v["next_event_in"] = timeOffset;
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

    nlohmann::json ObjectTracker::getUpcomingSatellitePassesWithPredictions()
    {
        nlohmann::json v;
        upcoming_satellite_passes_mtx.lock();
        v["upcoming_satellite_passes"] = upcoming_satellite_passes;
        v["upcoming_satellite_pass_points"] = upcoming_satellite_pass_points;

        std::vector<SatAzEl> upcoming_satellite_current_pos;
        double current_time = getTime();
        for(const auto& [norad, next_aos_time, next_los_time, max_elevation] : upcoming_satellite_passes)
        {
            SatAzEl sat_current_pos;
            if (norad_to_sat_object.count(norad))
            {
                // calculate only within time window
                if (current_time < next_los_time && current_time > next_aos_time)
                {
                    const predict_orbital_elements_t *satellite_object = norad_to_sat_object[norad]; // FIXME use stored sat objects!
                    predict_position satellite_orbit;
                    predict_observation satellite_observation_pos;
                    predict_orbit(satellite_object, &satellite_orbit, predict_to_julian_double(current_time));

                    general_mutex.lock(); // NOTE lock just in case observer station changes?
                    predict_observe_orbit(satellite_observer_station, &satellite_orbit, &satellite_observation_pos);
                    general_mutex.unlock();

                    sat_current_pos.az = satellite_observation_pos.azimuth * RAD_TO_DEG;
                    sat_current_pos.el = satellite_observation_pos.elevation * RAD_TO_DEG;
                }
            } else
            {
                logger->warn("Could not find norad in registry! norad:"+std::to_string(norad));
            }

            upcoming_satellite_current_pos.push_back(sat_current_pos);
        }

        v["upcoming_satellite_current_pos"] = upcoming_satellite_current_pos;
        upcoming_satellite_passes_mtx.unlock();
        return v;
    }
}
