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
    int hz;
    uint32_t iz;
    const float r = 3.442620;
    float value;
    float x;
    float y;

    hz = (int)shr3_seeded(jsr);
    iz = (hz & 127);

    if (fabs(hz) < kn[iz])
    {
        value = (float)(hz)*wn[iz];
    }
    else
    {
        for (; ; )
        {
            if (iz == 0)
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

                if (hz <= 0)
                {
                    value = -r - x;
                }
                else
                {
                    value = +r + x;
                }
                break;
            }

            x = (float)(hz)*wn[iz];

            if (fn[iz] + r4_uni(jsr) * (fn[iz - 1] - fn[iz])
                < exp(-0.5 * x * x))
            {
                value = x;
                break;
            }

            hz = (int)shr3_seeded(jsr);
            iz = (hz & 127);

            if (fabs(hz) < kn[iz])
            {
                value = (float)(hz)*wn[iz];
                break;
            }
        }
    }

    return value;
}

void r4_nor2_setup(uint32_t kn[128], float fn[128], float wn[128])
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

float r4_nor2(uint32_t* jsr, uint32_t kn[128], float fn[128], float wn[128], double* weightGrad = 0, double* biasGrad = 0, double* weight = 0, double* bias = 0)
{
    uint32_t temp;
    int32_t hz;
    uint8_t iz;
    float x, y;

    temp = *jsr;
    *jsr = (*jsr ^ (*jsr << 13));
    *jsr = (*jsr ^ (*jsr >> 17));
    *jsr = (*jsr ^ (*jsr << 5));
	hz = temp + *jsr;
    iz = hz & 127;
    if (hz < int32_t(kn[iz]) && ~hz + 1 < int32_t(kn[iz]))
        return wn[iz] * hz;
    
    for (;;)
    {
        if (iz == 0)
        {
            for (;;)
            {
                temp = *jsr;
                *jsr = (*jsr ^ (*jsr << 13));
                *jsr = (*jsr ^ (*jsr >> 17));
                *jsr = (*jsr ^ (*jsr << 5));
                x = (temp + *jsr) * 2.3283064365386963e-10f;
                x = (*(int32_t*)&x - 0x3f800000) * -2.3935259956e-8f;
                
                temp = *jsr;
                *jsr = (*jsr ^ (*jsr << 13));
                *jsr = (*jsr ^ (*jsr >> 17));
                *jsr = (*jsr ^ (*jsr << 5));
                y = (temp + *jsr) * 2.3283064365386963e-10f;
                y = (*(int32_t*)&y - 0x3f800000) * -8.24e-8f;
                
                if (x * x <= y + y)
                {
                    x += 3.442620f;
                    temp = hz & 0x80000000 ^ *(uint32_t*)&x;
                    return *(float*)&temp;
                }
            }
        }

        temp = *jsr;
        *jsr = (*jsr ^ (*jsr << 13));
        *jsr = (*jsr ^ (*jsr >> 17));
        *jsr = (*jsr ^ (*jsr << 5));
        y = (temp + *jsr) * 2.3283064365386963e-10f * (fn[iz - 1] - fn[iz]) + fn[iz];
        x = wn[iz] * hz;
        
        temp = 12338043.947595f * (-0.5f * x * x) + 1065101738.388696f;
		//printf("expected: %f, actual: %f\n", exp(-0.5f * x * x), *(float*)&temp);
        float res = exp(-0.5f * x * x);
		*biasGrad += (res - *(float*)&temp);
		*weightGrad += (res - *(float*)&temp) * (-0.5f * x * x);
        
        if (y < *(float*)&temp)
            return x;

        temp = *jsr;
        *jsr = (*jsr ^ (*jsr << 13));
        *jsr = (*jsr ^ (*jsr >> 17));
        *jsr = (*jsr ^ (*jsr << 5));
        hz = temp + *jsr;
        iz = hz & 127;
        if (hz < int32_t(kn[iz]) && ~hz + 1 < int32_t(kn[iz]))
            return wn[iz] * hz;
    }
}

float InvSqrt(float number)
{
    long i = 0x5F1FFFF9 - (*(long*)&number >> 1);
    float tmp = *(float*)&i;
    return tmp * 0.703952253f * (2.38924456f - number * tmp * tmp);
}

int main()
{
    /*uint32_t seed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    uint32_t tempSeed;
    float uniformRand;
    for (uint32_t i = 40; i--;)
    {
		tempSeed = seed1;
		seed1 = seed1 ^ (seed1 << 13);
		seed1 = seed1 ^ (seed1 >> 17);
		seed1 = seed1 ^ (seed1 << 5);
		uniformRand = (tempSeed + seed1) * 2.3283064365386963e-10f * 2 - 1;
        tempSeed = 12338084.563634f * uniformRand + 1065101642.533682f;
		printf("exp of %f is %f. an approximation is %f\n", uniformRand, exp(uniformRand), *(float*)&tempSeed);
    }
    return 0;*/
    
    uint32_t kn[128];
    float fn[128];
    float wn[128];
    
    uint32_t kn2[128];
    float fn2[128];
    float wn2[128];
    
    r4_nor_setup(kn, fn, wn);
	r4_nor2_setup(kn2, fn2, wn2);

    const uint32_t warmups = 20;
    const uint32_t loops = 10;
    const uint32_t bins = 128;
    const uint32_t samples = 100000;
    const float scale = 1000.0f / samples;
    const float min = -6.0f;
    const float max = 6.0f;
    const float bin_width = (max - min) / bins;
    
    /*uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
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
        printf("%f\n", fabs(scale * hist[i] - scale * hist2[i]));
        std::string spaces(scale * hist[i], ' ');
        printf("\t%s*%f\n", spaces.c_str(), scale * hist[i]);
		spaces = std::string(scale * hist2[i], ' ');
		printf("\t%s*%f\n", spaces.c_str(), scale * hist2[i]);
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
    printf("Time taken: %f microseconds\n", averageTime / loops);

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
    printf("Time taken: %f microseconds\n", averageTime / loops);*/

    double weightGrad, biasGrad;
    double weight = 12338043.947595;
	double bias = 1065101738.388696;
    for (;;)
    {
		weightGrad = 0;
		biasGrad = 0;
        uint32_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        for (uint32_t i = samples; i--;)
			r4_nor2(&seed, kn, fn, wn, &weightGrad, &biasGrad, &weight, &bias);
        /*printf("Weight Grad: %f\n", weightGrad);
		printf("Bias Grad: %f\n", biasGrad);*/
		weight += weightGrad * 0.01f;
        bias += biasGrad * 0.01f;
        printf("Weight: %f\n", weight);
        printf("Bias: %f\n", bias);
    }

    return 0;
}