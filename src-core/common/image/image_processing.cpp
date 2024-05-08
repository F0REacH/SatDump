#include "image_utils.h"
#include <algorithm>
#include <cmath>

namespace image
{
    int percentile(int *array, int size, float percentile)
    {
        float number_percent = (size + 1) * percentile / 100.0f;
        if (number_percent <= 1) // ==
            return array[0];
        else if (number_percent >= size) // ==
            return array[size - 1];
        else
            return array[(int)number_percent - 1] + (number_percent - (int)number_percent) * (array[(int)number_percent] - array[(int)number_percent - 1]);
    }

    void white_balance(Image &img, float percentileValue)
    {
        const float maxVal = img.maxval();
        const size_t d_height = img.height();
        const size_t d_width = img.width();

        int *sorted_array = new int[d_height * d_width];

        for (int c = 0; c < img.channels(); c++)
        {
            // Load the whole image band into our array
            for (size_t i = 0; i < d_height * d_width; i++)
                sorted_array[i] = img.get(c, i);

            // Sort it
            std::sort(&sorted_array[0], &sorted_array[d_width * d_height]);

            // Get percentiles
            int percentile1 = percentile(sorted_array, d_width * d_height, percentileValue);
            int percentile2 = percentile(sorted_array, d_width * d_height, 100.0f - percentileValue);

            for (size_t i = 0; i < d_width * d_height; i++)
            {
                long balanced = (img.get(c, i) - percentile1) * maxVal / (percentile2 - percentile1);
                if (balanced < 0)
                    balanced = 0;
                else if (balanced > maxVal)
                    balanced = maxVal;
                img.set(c, i, balanced);
            }
        }

        delete[] sorted_array;
    }

    void median_blur(Image &img)
    {
        for (int c = 0; c < img.channels(); c++)
        {

            int h = img.height();
            int w = img.width();

            std::vector<int> values(5);

            for (int x = 0; x < h; x++)
            {
                for (int y = 0; y < w; y++)
                {
                    values[0] =
                        values[1] =
                            values[2] =
                                values[3] =
                                    values[4] = img.get(c, x * w + y);

                    if (x != 0)
                        values[1] = img.get(c, (x - 1) * w + y);
                    if (y != 0)
                        values[2] = img.get(c, x * w + (y - 1));

                    if (x != h - 1)
                        values[3] = img.get(c, (x + 1) * w + y);
                    if (y != w - 1)
                        values[4] = img.get(c, x * w + (y + 1));

                    std::sort(values.begin(), values.end());

                    img.set(c, x * w + y, values[2]);
                }
            }
        }
    }

    void equalize(Image &img, bool per_channel)
    {
        for (int c = 0; c < (per_channel ? img.channels() : 1); c++)
        {
            if (c == 3) // Do not individual equalize alpha channel
                break;

            const int nlevels = img.maxval() + 1;
            size_t size = img.width() * img.height() * (per_channel ? 1 : img.channels());

            // Init histogram buffer
            int *histogram = new int[nlevels];
            for (int i = 0; i < nlevels; i++)
                histogram[i] = 0;

            // Compute histogram
            for (int px = 0; px < size; px++)
                histogram[img.get(c, px)]++;

            // Cummulative histogram
            int *cummulative_histogram = new int[nlevels];
            cummulative_histogram[0] = histogram[0];
            for (int i = 1; i < nlevels; i++)
                cummulative_histogram[i] = histogram[i] + cummulative_histogram[i - 1];

            // Scaling
            int *scaling = new int[nlevels];
            for (int i = 0; i < nlevels; i++)
                scaling[i] = round(cummulative_histogram[i] * (float(nlevels - 1) / size));

            // Apply
            for (int px = 0; px < size; px++)
                img.set(c, px, img.clamp(scaling[img.get(c, px)]));

            // Cleanup
            delete[] cummulative_histogram;
            delete[] scaling;
            delete[] histogram;
        }
    }

    void normalize(Image &img)
    {
        int max = 0;
        int min = img.maxval();

        // Get min and max
        for (size_t i = 0; i < img.size(); i++)
        {
            int val = img.get(i);

            if (val > max)
                max = val;
            if (val < min)
                min = val;
        }

        if (abs(max - min) == 0) // Avoid division by 0
            return;

        // Compute scaling factor
        float factor = img.maxval() / float(max - min);

        // Scale entire image
        for (size_t i = 0; i < img.size(); i++)
            img.set(i, img.clamp((img.get(i) - min) * factor));
    }

    void linear_invert(Image &img)
    {
        for (size_t i = 0; i < img.size(); i++)
            img.set(i, img.maxval() - img.get(i));
    }
}