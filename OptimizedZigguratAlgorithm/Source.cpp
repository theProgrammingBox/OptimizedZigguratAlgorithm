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

//uint32_t shr3_seeded(uint32_t* jsr)
//{
//    uint32_t value = *jsr;
//    *jsr = (*jsr ^ (*jsr << 13));
//    *jsr = (*jsr ^ (*jsr >> 17));
//    *jsr = (*jsr ^ (*jsr << 5));
//    return value + *jsr;
//}

float r4_uni(uint32_t* jsr)
{
    uint32_t jsr_input = *jsr;
    *jsr = (*jsr ^ (*jsr << 13));
    *jsr = (*jsr ^ (*jsr >> 17));
    *jsr = (*jsr ^ (*jsr << 5));
    return (jsr_input + *jsr) * 2.3283064365386963e-10f;
}

/*
TODO
0. see if kn can be int32_t
*/

float r4_nor(uint32_t* jsr, const uint32_t kn[128], const float fn[128], const float wn[128])
{
    int32_t randomInt32;
    uint32_t sevenBits;
    const float r = 3.442620f;
    float x, y;

    randomInt32 = *jsr;
    *jsr = (*jsr ^ (*jsr << 13));
    *jsr = (*jsr ^ (*jsr >> 17));
    *jsr = (*jsr ^ (*jsr << 5));
    randomInt32 += *jsr;
    sevenBits = randomInt32 & 127;
    if (randomInt32 < int32_t(kn[sevenBits]) && (~randomInt32 + 1) < int32_t(kn[sevenBits]))
        return wn[sevenBits] * randomInt32;
    /*if ((randomInt32 < int32_t(kn[sevenBits]) && (~randomInt32 + 1) < int32_t(kn[sevenBits])) != (fabs(randomInt32) < kn[sevenBits]))
        printf("Error, %d, %f, %d, %d, %d, %d\n", randomInt32, fabs(randomInt32), kn[sevenBits], (~randomInt32 + 1), int32_t(randomInt32) < kn[sevenBits], (~randomInt32 + 1) < kn[sevenBits]);*/
    
    for (;;)
    {
        if (sevenBits == 0)
        {
            for (;;)
            {
                x = -0.2904764 * log(r4_uni(jsr));
                y = -log(r4_uni(jsr));
                if (x * x <= y + y)
                {
                    if (randomInt32 <= 0)
                        return -r - x;
                    return +r + x;
                }
            }
        }

        x = wn[sevenBits] * randomInt32;
        if (fn[sevenBits] + r4_uni(jsr) * (fn[sevenBits - 1] - fn[sevenBits]) < exp(-0.5 * x * x))
            return x;

        randomInt32 = *jsr;
        *jsr = (*jsr ^ (*jsr << 13));
        *jsr = (*jsr ^ (*jsr >> 17));
        *jsr = (*jsr ^ (*jsr << 5));
        randomInt32 += *jsr;
        sevenBits = randomInt32 & 127;
        if (randomInt32 < int32_t(kn[sevenBits]) && (~randomInt32 + 1) < int32_t(kn[sevenBits]))
            return wn[sevenBits] * randomInt32;
    }
}

void randomNormalSetup(uint32_t kn[128], float fn[128], float wn[128])
{
    double dn = 3.442619855899;
    const double m1 = 2147483648.0;
    const double vn = 9.91256303526217e-03;
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

float randomNormal(uint32_t& seed, const uint32_t kn[128], const float fn[128], const float wn[128])
{
    uint32_t tempSeed;
    int32_t randomInt;
    uint8_t sevenBits;
    float x, y;

    for (;;)
    {
        tempSeed = seed;
        seed = (seed ^ (seed << 13));
        seed = (seed ^ (seed >> 17));
        seed = (seed ^ (seed << 5));
        randomInt = tempSeed + seed;
        sevenBits = randomInt & 0x7f;
        if (randomInt < kn[sevenBits] || randomInt & 0x80000000 && ~randomInt + 1 < kn[sevenBits])
            return wn[sevenBits] * randomInt;

        if (sevenBits == 0)
        {
            for (;;)
            {
                tempSeed = seed;
                seed = (seed ^ (seed << 13));
                seed = (seed ^ (seed >> 17));
                seed = (seed ^ (seed << 5));
                x = (tempSeed + seed) * 2.32830643654e-10f;
                x = (*(int32_t*)&x - 0x3f800000) * -2.3935259956e-8f;

                tempSeed = seed;
                seed = (seed ^ (seed << 13));
                seed = (seed ^ (seed >> 17));
                seed = (seed ^ (seed << 5));
                y = (tempSeed + seed) * 2.32830643654e-10f;
                y = (*(int32_t*)&y - 0x3f800000) * -8.24e-8f;

                if (x * x <= y + y)
                {
                    x += 3.442620f;
                    tempSeed = *(uint32_t*)&x ^ randomInt & 0x80000000;
                    return *(float*)&tempSeed;
                }
            }
        }

        x = wn[sevenBits] * randomInt;
        tempSeed = seed;
        seed = (seed ^ (seed << 13));
        seed = (seed ^ (seed >> 17));
        seed = (seed ^ (seed << 5));
        if (fn[sevenBits] + (tempSeed + seed) * 2.32830643654e-10f * (fn[sevenBits - 1] - fn[sevenBits]) < exp(-0.5f * x * x))
            return x;
    }
}

int main()
{
    uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    uint32_t kn[128];
    float fn[128];
    float wn[128];

    //randomNormalSetup(kn, fn, wn);
    r4_nor_setup(kn, fn, wn);

    const uint32_t bins = 54;
    const uint32_t samples = 100000000;
    const float scale = 1000.0f / samples;
    const float min = -3.0f;
    const float max = 3.0f;
    const float bin_width = (max - min) / bins;
    
	auto start = std::chrono::high_resolution_clock::now();
    uint32_t hist[bins];
    memset(hist, 0, sizeof(hist));
    for (uint32_t i = 0; i < samples; i++)
    {
        //uint32_t bin = (randomNormal(seed, kn, fn, wn) - min) / bin_width;
        uint32_t bin = (r4_nor(&seed, kn, fn, wn) - min) / bin_width;
        if (bin < bins && bin >= 0)
            hist[bin]++;
    }
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	printf("Time taken: %lld microseconds\n", duration.count());

    printf("----------------------------------------\n");
    for (uint32_t i = 0; i < bins; i++)
    {
        for (uint32_t j = scale * hist[i]; j--;)
            printf("*");
        printf("\n");
    }
    printf("----------------------------------------\n");

    return 0;
}