Texture2D src : register(t0);
SamplerState s;
cbuffer constant0 : register(b0) {
    float2 resolution;
    float2 originalCenter;
    float amplitude;
    float inverseRadius;
    float phase;
    float shape;
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

float4 psmain(float4 pos : SV_Position) : SV_Target {
    float2 p = pos.xy - originalCenter;
    float wave = calculateWave(length(p) * inverseRadius + phase, (int)(shape + 0.5));
    float sine, cosine;
    sincos(wave * amplitude, sine, cosine);
    float2 sourcePosition = float2(cosine * p.x + sine * p.y, -sine * p.x + cosine * p.y);
    float2 uv = (originalCenter + sourcePosition) / resolution;
    if (any(uv < 0.0) || any(uv > 1.0)) return float4(0, 0, 0, 0);
    return src.Sample(s, uv);
}
