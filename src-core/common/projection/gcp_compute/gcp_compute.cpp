#include "gcp_compute.h"
#include "nlohmann/json_utils.h"
#include "../sat_proj/sat_proj.h"

#include "logger.h"

namespace satdump
{
    namespace gcp_compute
    {
        std::vector<satdump::projection::GCP> compute_gcps(nlohmann::ordered_json cfg, TLE tle, nlohmann::ordered_json timestamps, int width, int height)
        {
            std::vector<satdump::projection::GCP> gcps;

            std::shared_ptr<SatelliteProjection> projection = get_sat_proj(cfg, tle, timestamps);

            std::vector<int> values;
            for (int x = 0; x < projection->img_size_x; x += projection->gcp_spacing_x)
                values.push_back(x);
            values.push_back(projection->img_size_x - 1);

            geodetic::geodetic_coords_t position;
            double ratio_x = 1, ratio_y = 1;

            if (width != -1)
                ratio_x = round((double)projection->img_size_x / (double)width);
            if (height != -1)
                ratio_y = round((double)projection->img_size_y / (double)height);

            logger->critical("{:d}x{:d}  {:d}x{:d}    {:f} {:f}", projection->img_size_x, projection->img_size_y, width, height, ratio_x, ratio_y);

            bool last_was_invalid = false;
            for (int y = 0; y < projection->img_size_y; y++)
            {
                for (int x : values)
                {
                    if (y % projection->gcp_spacing_y == 0 || y + 1 == (int)timestamps.size() || last_was_invalid)
                    {
                        if (projection->get_position(x * ratio_x, y * ratio_y, position))
                        {
                            last_was_invalid = true;
                            continue;
                        }

                        gcps.push_back({(double)x, (double)y, (double)position.lon, (double)position.lat});
                    }

                    last_was_invalid = false;
                }
            }

            return gcps;
        }
    }
}