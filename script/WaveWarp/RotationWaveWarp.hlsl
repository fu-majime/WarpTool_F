Texture2D src : register(t0);
SamplerState s;
cbuffer constant0 : register(b0) {
    float2 resolution;
    float2 originalCenter;
    float amplitude;
    float inverseRadius;
    float phase;
    float shape;
    float reference;
    float randHeight;
    float randWidth;
    float randSeed;
    float randTime;
};
static const float TAU = 6.28318530718;

static const int WAVE_SINE = 1;
static const int WAVE_SQUARE = 2;
static const int WAVE_TRIANGLE = 3;
static const int WAVE_SAWTOOTH_POS = 4;
static const int WAVE_SAWTOOTH_NEG = 5;
static const int WAVE_CIRCLE = 6;
static const int WAVE_SEMICIRCLE_POS = 7;
static const int WAVE_SEMICIRCLE_NEG = 8;

float calculateWave(float input, int waveType) {
    switch (waveType) {
        case WAVE_SINE:
            return sin(TAU * input);
        case WAVE_SQUARE:
            return sign(sin(TAU * input));
        case WAVE_TRIANGLE:
            return (1.0 - 4.0 * abs(0.5 - frac(input + 0.25)));
        case WAVE_SAWTOOTH_POS:
            return 2.0 * frac(input) - 1.0;
        case WAVE_SAWTOOTH_NEG:
            return 1.0 - 2.0 * frac(input);
        case WAVE_CIRCLE: {
            float r = frac(input);
            float val1 = -16.0 * r * r + 8.0 * r;
            float val2 = -16.0 * r * r + 24.0 * r - 8.0;
            if (val1 > 0) {
                return sqrt(val1);
            } else if (val2 > 0) {
                return -sqrt(val2);
            } else {
                return 0;
            }
        }
        case WAVE_SEMICIRCLE_POS: {
            float r = frac(input);
            return sqrt(-16.0 * r * r + 16.0 * r) - 1;
        }
        case WAVE_SEMICIRCLE_NEG: {
            float r = frac(input);
            return 1 - sqrt(-16.0 * r * r + 16.0 * r);
        }
        default:
            return sin(TAU * input);
    }
}

float animatedHash(int n, int a, int b, int c) {
    int step = (int)floor(randTime);
    int2 v = n + (int)randSeed + int2(step, step + 1) * 104729;
    v = (v << 13) ^ v;
    float2 values = (float2)((v * (v * v * a + b) + c) & 0x7fffffff) / 2147483647.0;
    float t = frac(randTime);
    return lerp(values.x, values.y, t * t * (3.0 - 2.0 * t));
}
float hash1(int n) { return animatedHash(n, 17389, 611953, 1611623773); }
float hash2(int n) { return animatedHash(n, 27449, 746773, 1824261409); }

// ランダム値の累積数列の近似
float stepNoise(int n, float amplitude) {
    float t = (float)n / 2.0 + hash2(n) / 2.0;
    return lerp((float)n, t, amplitude);
}

// stepNoise(index) <= x < stepNoise(index+1) を満たすindexを見つける
int estimateStepNoiseIndex(float x, float amplitude) {
    int indexApprox = (int)floor(x / (1.0 - amplitude * 0.5));
    float stepNoiseValue = stepNoise(indexApprox, amplitude);
    int isGreater = (stepNoiseValue > x) ? 1 : 0;
    return indexApprox - isGreater;
}

float applyWaveWidthRandomness(float input, float amount, int waveType) {
    switch (waveType) {
        case WAVE_SQUARE:
        case WAVE_SEMICIRCLE_POS:
        case WAVE_SEMICIRCLE_NEG: {
            float scaledInput = 2.0 * input;
            int index = estimateStepNoiseIndex(scaledInput, amount);
            return 0.5 * index
            + 0.5 * (scaledInput - stepNoise(index, amount))
            / (stepNoise(index + 1, amount) - stepNoise(index, amount));
        }
        case WAVE_SAWTOOTH_POS:
        case WAVE_SAWTOOTH_NEG: {
            int index = estimateStepNoiseIndex(input, amount);
            return index
            + (input - stepNoise(index, amount))
            / (stepNoise(index + 1, amount) - stepNoise(index, amount));
        }
        case WAVE_SINE:
        case WAVE_TRIANGLE:
        case WAVE_CIRCLE: {
            float scaledInput = 2.0 * input + 0.5;
            int index = estimateStepNoiseIndex(scaledInput, amount);
            return 0.5 * index - 0.25
            + 0.5 * (scaledInput - stepNoise(index, amount))
            / (stepNoise(index + 1, amount) - stepNoise(index, amount));
        }
        default:
            return input;
    }
}

float applyWaveHeightRandomness(float wave, float input, float amount, int waveType) {
    switch (waveType) {
        case WAVE_SQUARE: {
            float amp = 1.0 - hash1((int)floor(2 * input)) * amount;
            return wave * amp;
        }
        case WAVE_SAWTOOTH_POS:
        case WAVE_SAWTOOTH_NEG: {
            float amp = 1.0 - hash1(2 * (int)floor(input)) * amount;
            return wave * amp;
        }
        case WAVE_SINE:
        case WAVE_TRIANGLE:
        case WAVE_CIRCLE: {
            float amp0 = 0, amp1 = 0;
            if (frac(input + 0.25) < 0.5) {
                amp0 = 1.0 - hash1((int)floor(2 * input - 0.5)) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input + 0.5)) * amount;
                return ((amp1 + amp0) * wave + (amp1 - amp0)) * 0.5;
            } else {
                amp0 = 1.0 - hash1((int)floor(2 * input - 0.5)) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input + 0.5)) * amount;
                return ((amp1 + amp0) * wave + (amp0 - amp1)) * 0.5;
            }
        }
        case WAVE_SEMICIRCLE_POS:
        case WAVE_SEMICIRCLE_NEG: {
            float amp0 = 0, amp1 = 0;
            if (frac(input) < 0.5) {
                amp0 = 1.0 - hash1((int)floor(2 * input - 1.0)) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input)) * amount;
                return ((amp1 + amp0) * wave + (amp1 - amp0)) * 0.5;
            } else {
                amp0 = 1.0 - hash1((int)floor(2 * input - 1.0)) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input)) * amount;
                return ((amp1 + amp0) * wave + (amp0 - amp1)) * 0.5;
            }
        }
        default:
            return wave;
    }
}

float4 psmain(float4 pos : SV_Position) : SV_Target {
    float2 p = pos.xy - originalCenter;
    int waveType = (int)(shape + 0.5);
    float input = applyWaveWidthRandomness(length(p) * inverseRadius - phase, randWidth * 0.01, waveType);
    float wave = calculateWave(input, waveType);
    wave = applyWaveHeightRandomness(wave, input, randHeight * 0.01, waveType);
    float sine, cosine;
    sincos((wave + reference) * amplitude, sine, cosine);
    float2 sourcePosition = float2(cosine * p.x - sine * p.y, sine * p.x + cosine * p.y);
    float2 uv = (originalCenter + sourcePosition) / resolution;
    if (any(uv < 0.0) || any(uv > 1.0)) return float4(0, 0, 0, 0);
    return src.Sample(s, uv);
}
