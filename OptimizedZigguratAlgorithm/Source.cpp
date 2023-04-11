#include <iostream>
#include <chrono>

void r4_nor_setup(uint32_t kn[128], float fn[128], float wn[128])
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
    int32_t hz;
    uint32_t iz;
    float x;
    float y;

    hz = (int)shr3_seeded(jsr);
    iz = (hz & 127);

    if (fabs(hz) < kn[iz])
        return (float)(hz)*wn[iz];
    
    for (;;)
    {
        if (iz == 0)
        {
            for (;;)
            {
                x = -0.2904764 * log(r4_uni(jsr));
                y = -log(r4_uni(jsr));
                if (x * x <= y + y)
                {
                    if (hz <= 0)
                        return -3.442620f - x;
                    return +3.442620f + x;
                }
            }
        }

        x = (float)(hz)*wn[iz];

        if (fn[iz] + r4_uni(jsr) * (fn[iz - 1] - fn[iz]) < exp(-0.5 * x * x))
            return x;

        hz = (int)shr3_seeded(jsr);
        iz = (hz & 127);

        if (fabs(hz) < kn[iz])
            return (float)(hz)*wn[iz];
    }
}

/*
//TODO
//0. see if kn can be int32_t
//1. check wiki optimizations
//2. unroll loop, get stats on how many iterations it takes

float r4_nor(uint32_t* seed, const int32_t kn[128], const float fn[128], const float wn[128])
{
    int32_t randomInt32;
    uint8_t sevenBits;
    float x, y;

    randomInt32 = *seed;
    *seed = (*seed ^ (*seed << 13));
    *seed = (*seed ^ (*seed >> 17));
    *seed = (*seed ^ (*seed << 5));
    randomInt32 += *seed;
    sevenBits = randomInt32 & 127;
    if (randomInt32 < kn[sevenBits] && ~randomInt32 + 1 < kn[sevenBits])
        return wn[sevenBits] * randomInt32;
    //if ((randomInt32 < int32_t(kn[sevenBits]) && (~randomInt32 + 1) < int32_t(kn[sevenBits])) != (fabs(randomInt32) < kn[sevenBits]))
        //printf("Error, %d, %f, %d, %d, %d, %d\n", randomInt32, fabs(randomInt32), kn[sevenBits], (~randomInt32 + 1), int32_t(randomInt32) < kn[sevenBits], (~randomInt32 + 1) < kn[sevenBits]);

    for (;;)
    {
        if (sevenBits == 0)
        {
            for (;;)
            {
				// randomInt32 is an int, so range is kinda messed up
                randomInt32 = *seed;
                *seed = (*seed ^ (*seed << 13));
                *seed = (*seed ^ (*seed >> 17));
                *seed = (*seed ^ (*seed << 5));
				x = uint32_t(randomInt32 + *seed) * 2.3283064365386963e-10f;
                //x = (*(int32_t*)&x - 0x3f800000) * -2.3935259956e-8f;
                x = -0.2904764f * log(x);

                randomInt32 = *seed;
				*seed = (*seed ^ (*seed << 13));
				*seed = (*seed ^ (*seed >> 17));
				*seed = (*seed ^ (*seed << 5));
				y = uint32_t(randomInt32 + *seed) * 2.3283064365386963e-10f;
                //y = (*(int32_t*)&y - 0x3f800000) * -8.24e-8f;
				y = -log(y);
                
                if (x * x <= y + y)
                {
                    x += 3.442620f;
                    randomInt32 = *(uint32_t*)&x ^ randomInt32 & 0x80000000;
                    return *(float*)&randomInt32;
                }
            }
        }

        randomInt32 = *seed;
        *seed = (*seed ^ (*seed << 13));
        *seed = (*seed ^ (*seed >> 17));
        *seed = (*seed ^ (*seed << 5));
        y = (randomInt32 + *seed) * 2.3283064365386963e-10f;
        x = wn[sevenBits] * randomInt32;
        if (fn[sevenBits] + y * (fn[sevenBits - 1] - fn[sevenBits]) < exp(-0.5f * x * x))
            return x;

        randomInt32 = *seed;
        *seed = (*seed ^ (*seed << 13));
        *seed = (*seed ^ (*seed >> 17));
        *seed = (*seed ^ (*seed << 5));
        randomInt32 += *seed;
        sevenBits = randomInt32 & 127;
        if (randomInt32 < kn[sevenBits] && ~randomInt32 + 1 < kn[sevenBits])
            return wn[sevenBits] * randomInt32;
    }
}*/

int main()
{
    uint32_t kn[128];
    float fn[128];
    float wn[128];
    
    r4_nor_setup(kn, fn, wn);

    const uint32_t bins = 54;
    const uint32_t samples = 40000000;
    const float scale = 1000.0f / samples;
    const float min = -3.0f;
    const float max = 3.0f;
    const float bin_width = (max - min) / bins;

	float averageTime = 0.0f;
    for (uint32_t j = 10; j--;)
    {
        uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        auto start = std::chrono::high_resolution_clock::now();
        for (uint32_t i = samples; i--;)
            r4_nor(&seed, kn, fn, wn);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        averageTime += duration.count() * 1e-3f;
    }
	printf("Time taken: %f microseconds\n", averageTime * 0.1f);

    uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    uint32_t hist[bins];
    memset(hist, 0, sizeof(hist));
    for (uint32_t i = 0; i < samples; i++)
    {
        //uint32_t bin = (randomNormal(seed, kn, fn, wn) - min) / bin_width;
        uint32_t bin = (r4_nor(&seed, kn, fn, wn) - min) / bin_width;
        if (bin < bins && bin >= 0)
            hist[bin]++;
    }

    printf("----------------------------------------\n");
    for (uint32_t i = 0; i < bins; i++)
    {
        /*for (uint32_t j = scale * hist[i]; j--;)
            printf("*");
        printf("\n");*/
        
		std::string spaces(scale * hist[i], ' ');
		printf("%s%f\n", spaces.c_str(), scale * hist[i]);
    }
    printf("----------------------------------------\n");

    return 0;
}