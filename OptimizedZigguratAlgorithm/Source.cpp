#include <iostream>
#include <chrono>

void r4_nor_setup(uint32_t kn[128], float fn[128], float wn[128])
{
    double dn = 3.442619855899;
    int i;
    const double m1 = 2147483648.0;
    double q;
    double tn = 3.442619855899;
    const double vn = 9.91256303526217E-03;

    q = vn / exp(-0.5 * dn * dn);

    kn[0] = (uint32_t)((dn / q) * m1);
    kn[1] = 0;

    wn[0] = (float)(q / m1);
    wn[127] = (float)(dn / m1);

    fn[0] = 1.0;
    fn[127] = (float)(exp(-0.5 * dn * dn));

    for (i = 126; 1 <= i; i--)
    {
        dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
        kn[i + 1] = (uint32_t)((dn / tn) * m1);
        tn = dn;
        fn[i] = (float)(exp(-0.5 * dn * dn));
        wn[i] = (float)(dn / m1);
    }

    return;
}

uint32_t shr3_seeded(uint32_t* jsr)
{
    uint32_t value;

    value = *jsr;

    *jsr = (*jsr ^ (*jsr << 13));
    *jsr = (*jsr ^ (*jsr >> 17));
    *jsr = (*jsr ^ (*jsr << 5));

    value = value + *jsr;

    return value;
}

float r4_uni(uint32_t* jsr)
{
    uint32_t jsr_input;
    float value;

    jsr_input = *jsr;

    *jsr = (*jsr ^ (*jsr << 13));
    *jsr = (*jsr ^ (*jsr >> 17));
    *jsr = (*jsr ^ (*jsr << 5));

    value = fmod(0.5
        + (float)(jsr_input + *jsr) / 65536.0 / 65536.0, 1.0);

    return value;
}

float r4_nor(uint32_t* jsr, uint32_t kn[128], float fn[128], float wn[128])
{
    int int32Seed;
    uint32_t int8Seed;
    const float r = 3.442620;
    float value;
    float x;
    float y;

    int32Seed = (int)shr3_seeded(jsr);
    int8Seed = (int32Seed & 127);

    if (fabs(int32Seed) < kn[int8Seed])
    {
        value = (float)(int32Seed)*wn[int8Seed];
    }
    else
    {
        for (; ; )
        {
            if (int8Seed == 0)
            {
                for (; ; )
                {
                    x = -0.2904764 * log(r4_uni(jsr));
                    y = -log(r4_uni(jsr));
                    if (x * x <= y + y)
                    {
                        break;
                    }
                }

                if (int32Seed <= 0)
                {
                    value = -r - x;
                }
                else
                {
                    value = +r + x;
                }
                break;
            }

            x = (float)(int32Seed)*wn[int8Seed];

            if (fn[int8Seed] + r4_uni(jsr) * (fn[int8Seed - 1] - fn[int8Seed])
                < exp(-0.5 * x * x))
            {
                value = x;
                break;
            }

            int32Seed = (int)shr3_seeded(jsr);
            int8Seed = (int32Seed & 127);

            if (fabs(int32Seed) < kn[int8Seed])
            {
                value = (float)(int32Seed)*wn[int8Seed];
                break;
            }
        }
    }

    return value;
}

void r4_nor2_setup(int32_t kn[128], float fn[128], float wn[128])
{
    double dn = 3.442619855899;
    const double m1 = 2147483648.0;
    const double vn = 9.91256303526217E-03;
    double q = vn / exp(-0.5 * dn * dn);

    kn[0] = dn / q * m1;
    kn[1] = 0;
    wn[0] = q / m1;
    wn[127] = dn / m1;
    fn[0] = 1.0;
    fn[127] = exp(-0.5 * dn * dn);

    double tn;
    for (uint8_t i = 126; 1 <= i; i--)
    {
        tn = dn;
        dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
        kn[i + 1] = dn / tn * m1;
        fn[i] = exp(-0.5 * dn * dn);
        wn[i] = dn / m1;
    }

    return;
}

float r4_nor2(uint32_t* jsr, int32_t kn[128], float fn[128], float wn[128])
{
    uint32_t uint32Temp;
    int32_t int32Seed;
    uint8_t int8Seed;
    float x, y;

    uint32Temp = *jsr;
    *jsr ^= *jsr << 13;
    *jsr ^= *jsr >> 17;
    *jsr ^= *jsr << 5;
	int32Seed = uint32Temp + *jsr;
    int8Seed = int32Seed & 127;
    if ((int32Seed ^ (int32Seed >> 31)) - (int32Seed >> 31) < kn[int8Seed])
        return wn[int8Seed] * int32Seed;
    
    for (;;)
    {
        if (int8Seed == 0)
        {
            for (;;)
            {
                uint32Temp = *jsr;
                *jsr ^= *jsr << 13;
                *jsr ^= *jsr >> 17;
                *jsr ^= *jsr << 5;
                x = (uint32Temp + *jsr) * 2.3283064365386963e-10f;
                x = (1065353216 - *(int32_t*)&x) * 2.30830217163e-08f;
                
                uint32Temp = *jsr;
                *jsr ^= *jsr << 13;
                *jsr ^= *jsr >> 17;
                *jsr ^= *jsr << 5;
                y = (uint32Temp + *jsr) * 2.3283064365386963e-10f;
                y = (1065353216 - *(int32_t*)&y) * 7.94660834913e-08f;
                
                if (x * x <= y + y)
                {
                    x += 3.442620f;
                    uint32Temp = int32Seed & 0x80000000 ^ *(uint32_t*)&x;
                    return *(float*)&uint32Temp;
                }
            }
        }

        uint32Temp = *jsr;
        *jsr ^= *jsr << 13;
        *jsr ^= *jsr >> 17;
        *jsr ^= *jsr << 5;
        y = (uint32Temp + *jsr) * 2.3283064365386963e-10f;
        x = wn[int8Seed] * int32Seed;
        uint32Temp = -6169045.423972f * x * x + 1065101626.864132f;
        if (y * (fn[int8Seed - 1] - fn[int8Seed]) + fn[int8Seed] < *(float*)&uint32Temp)
            return x;

        uint32Temp = *jsr;
        *jsr ^= *jsr << 13;
        *jsr ^= *jsr >> 17;
        *jsr ^= *jsr << 5;
        int32Seed = uint32Temp + *jsr;
        int8Seed = int32Seed & 127;
        if ((int32Seed ^ (int32Seed >> 31)) - (int32Seed >> 31) < kn[int8Seed])
            return wn[int8Seed] * int32Seed;
    }
}


int main()
{
    uint32_t kn[128];
    float fn[128];
    float wn[128];

    int32_t kn2[128];
    float fn2[128];
    float wn2[128];

    r4_nor_setup(kn, fn, wn);
    r4_nor2_setup(kn2, fn2, wn2);

    const uint32_t warmups = 20;
    const uint32_t loops = 20;
    const uint32_t bins = 128;
    const uint32_t samples = 10000000;
    const float scale = 1000.0f / samples;
    const float min = -6.0f;
    const float max = 6.0f;
    const float bin_width = (max - min) / bins;
    
    uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    uint32_t hist[bins];
    memset(hist, 0, sizeof(hist));
    for (uint32_t i = 0; i < samples; i++)
    {
        uint32_t bin = (r4_nor(&seed, kn, fn, wn) - min) / bin_width;
        if (bin < bins && bin >= 0)
            hist[bin]++;
    }

    seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    uint32_t hist2[bins];
    memset(hist2, 0, sizeof(hist2));
    for (uint32_t i = 0; i < samples; i++)
    {
        uint32_t bin = (r4_nor2(&seed, kn2, fn2, wn2) - min) / bin_width;
        if (bin < bins && bin >= 0)
            hist2[bin]++;
    }

    for (uint32_t i = 0; i < bins; i++)
    {
        /*printf("%f\%%%\n", (1.0f - float(hist2[i]) / hist[i]) * 100.0f);
        std::string spaces(scale * hist[i], ' ');
        printf("\t%s*%f\n", spaces.c_str(), scale * hist[i]);
        spaces = std::string(scale * hist2[i], ' ');
        printf("\t%s*%f\n", spaces.c_str(), scale * hist2[i]);*/
        std::string spaces(scale * hist2[i], ' ');
		printf("%s*\n", spaces.c_str());
    }

    for (uint32_t j = warmups; j--;)
    {
        uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        for (uint32_t i = samples; i--;)
            r4_nor(&seed, kn, fn, wn);
    }

    float averageTime = 0.0f;
    for (uint32_t j = loops; j--;)
    {
        uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        auto start = std::chrono::high_resolution_clock::now();
        for (uint32_t i = samples; i--;)
            r4_nor(&seed, kn, fn, wn);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        averageTime += duration.count() * 1e-3f;
    }
	double time1 = averageTime / loops;
	printf("Time taken: %f microseconds\n", time1);

    averageTime = 0.0f;
    for (uint32_t j = loops; j--;)
    {
        uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        auto start = std::chrono::high_resolution_clock::now();
        for (uint32_t i = samples; i--;)
            r4_nor2(&seed, kn2, fn2, wn2);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        averageTime += duration.count() * 1e-3f;
    }
	double time2 = averageTime / loops;
	printf("Time taken: %f microseconds\n", time2);

    printf("Speedup: %f%%%\n", (1 - time2 / time1) * 100);

    return 0;
}