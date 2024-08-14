#include "clock_recovery_mm.h"
#include "common/dsp/buffer.h"
#include "common/dsp/window/window.h"
#include "common/dsp/block.h"

#define DO_BRANCH 0

namespace ndsp
{
    template <typename T>
    ClockRecoveryMM<T>::ClockRecoveryMM()
        : ndsp::Block(std::is_same_v<T, float> ? "clock_recovery_mm_ff" : "clock_recovery_mm_cc", {{sizeof(T)}}, {{sizeof(T)}}),
          p_2T(0),
          p_1T(0),
          p_0T(0),
          c_2T(0),
          c_1T(0),
          c_0T(0)
    {
    }

    template <typename T>
    ClockRecoveryMM<T>::~ClockRecoveryMM()
    {
        if (buffer != nullptr)
            volk_free(buffer);
    }

    template <typename T>
    void ClockRecoveryMM<T>::start()
    {
        set_params();
        ndsp::buf::init_nafifo_stdbuf<T>(outputs[0], 2, ((ndsp::buf::StdBuf<T> *)inputs[0]->read_buf())->max); // TODO FIX
        ndsp::Block::start();
    }

    template <typename T>
    void ClockRecoveryMM<T>::set_params(nlohmann::json p)
    {
        if (p.contains("nfilt"))
            d_nfilt = p["nfilt"];
        if (p.contains("ntaps"))
            d_ntaps = p["ntaps"];
        if (p.contains("mu"))
            d_mu = p["mu"];
        if (p.contains("omega"))
            d_omega = p["omega"];
        if (p.contains("omega_gain"))
            d_omega_gain = p["omega_gain"];
        if (p.contains("mu_gain"))
            d_mu_gain = p["mu_gain"];
        if (p.contains("omega_relative_limit"))
            d_omega_relative_limit = p["omega_relative_limit"];

        // General setup
        mu = d_mu;
        omega = d_omega;
        omega_gain = d_omega_gain;
        mu_gain = d_mu_gain;
        omega_relative_limit = d_omega_relative_limit;

        // Omega setup
        omega_mid = omega;
        omega_limit = omega_relative_limit * omega;

        // Init interpolator
        pfb.init(dsp::windowed_sinc(d_nfilt * d_ntaps, dsp::hz_to_rad(0.5 / (double)d_nfilt, 1.0), dsp::window::nuttall, d_nfilt), d_nfilt);

        // Buffer
#if DO_BRANCH
        buffer = create_volk_buffer<T>(pfb.ntaps * 4);
#else
        buffer = dsp::create_volk_buffer<T>(((ndsp::buf::StdBuf<T> *)inputs[0]->read_buf())->max * 2);
#endif
    }

    template <typename T>
    void ClockRecoveryMM<T>::work()
    {
        if (!inputs[0]->read())
        {
            auto *rbuf = (ndsp::buf::StdBuf<T> *)inputs[0]->read_buf();
            auto *wbuf = (ndsp::buf::StdBuf<T> *)outputs[0]->write_buf();

            int nsamples = rbuf->cnt;

            // Copy NTAPS samples in the buffer from input, as that's required for the last samples
#if DO_BRANCH
            memcpy(&buffer[pfb.ntaps - 1], rbuf->dat, (pfb.ntaps - 1) * sizeof(T));
#else
            memcpy(&buffer[pfb.ntaps - 1], rbuf->dat, nsamples * sizeof(T));
#endif
            ouc = 0;

            for (; inc < nsamples && ouc < rbuf->max;)
            {
                if constexpr (std::is_same_v<T, complex_t>)
                {
                    // Propagate delay
                    p_2T = p_1T;
                    p_1T = p_0T;
                    c_2T = c_1T;
                    c_1T = c_0T;
                }

                // Compute output
                int imu = (int)rint(mu * pfb.nfilt);
                if (imu < 0) // If we're out of bounds, clamp
                    imu = 0;
                if (imu >= pfb.nfilt)
                    imu = pfb.nfilt - 1;

                if constexpr (std::is_same_v<T, float>)
                {
#if DO_BRANCH
                    if (inc < (pfb.ntaps - 1))
                        volk_32f_x2_dot_prod_32f(&sample, &buffer[inc], pfb.taps[imu], pfb.ntaps);
                    else
                        volk_32f_x2_dot_prod_32f(&sample, &rbuf->dat[inc - (pfb.ntaps - 1)], pfb.taps[imu], pfb.ntaps);
#else
                    volk_32f_x2_dot_prod_32f(&sample, &buffer[inc], pfb.taps[imu], pfb.ntaps);
#endif

                    // Phase error
                    phase_error = (last_sample < 0 ? -1.0f : 1.0f) * sample - (sample < 0 ? -1.0f : 1.0f) * last_sample;
                    phase_error = BRANCHLESS_CLIP(phase_error, 1.0);
                    last_sample = sample;

                    // Write output sample
                    wbuf->dat[ouc++] = sample;
                }
                if constexpr (std::is_same_v<T, complex_t>)
                {
#if DO_BRANCH
                    if (inc < (pfb.ntaps - 1))
                        volk_32fc_32f_dot_prod_32fc((lv_32fc_t *)&p_0T, (lv_32fc_t *)&buffer[inc], pfb.taps[imu], pfb.ntaps);
                    else
                        volk_32fc_32f_dot_prod_32fc((lv_32fc_t *)&p_0T, (lv_32fc_t *)&rbuf->dat[inc - (pfb.ntaps - 1)], pfb.taps[imu], pfb.ntaps);
#else
                    volk_32fc_32f_dot_prod_32fc((lv_32fc_t *)&p_0T, (lv_32fc_t *)&buffer[inc], pfb.taps[imu], pfb.ntaps);
#endif

                    // Slice it
                    c_0T = complex_t(p_0T.real > 0.0f ? 1.0f : 0.0f, p_0T.imag > 0.0f ? 1.0f : 0.0f);

                    // Phase error
                    phase_error = (((p_0T - p_2T) * c_1T.conj()) - ((c_0T - c_2T) * p_1T.conj())).real;
                    phase_error = BRANCHLESS_CLIP(phase_error, 1.0);

                    // Write output
                    wbuf->dat[ouc++] = p_0T;
                }

                // Adjust omega
                omega = omega + omega_gain * phase_error;
                omega = omega_mid + BRANCHLESS_CLIP((omega - omega_mid), omega_limit);

                // Adjust phase
                mu = mu + omega + mu_gain * phase_error;
                inc += int(floor(mu));
                mu -= floor(mu);
                if (inc < 0)
                    inc = 0;
            }

            inc -= nsamples;

            if (inc < 0)
                inc = 0;

                // We need some history for the next run, so copy it over into our buffer
#if DO_BRANCH
            memcpy(buffer, &rbuf->dat[nsamples - pfb.ntaps + 1], (pfb.ntaps - 1) * sizeof(T));
#else
            memmove(&buffer[0], &buffer[nsamples], pfb.ntaps * sizeof(T));
#endif

            wbuf->cnt = ouc;
            inputs[0]->flush();
            outputs[0]->write();
        }
    }

    template class ClockRecoveryMM<complex_t>;
    template class ClockRecoveryMM<float>;
}
