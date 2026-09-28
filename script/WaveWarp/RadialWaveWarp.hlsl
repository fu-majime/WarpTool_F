Texture2D src : register(t0);
SamplerState s;
cbuffer constant0 : register(b0) {
    float2 resolution;
    float2 originalCenter;
    float baseRadius;
    float amplitude;
    float frequency;
    float phase;
    float shape;
    float randHeight;
    float randWidth;
    float randSeed;
    float basePosition;
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

float animatedHash(int n, int a, int b, int c, float period) {
    float x = period > 0.0 ? frac(n / max(period, 1e-6)) * 64.0 : 0.0;
    int index = (int)floor(x);
    int2 indices = period > 0.0 ? int2(index, (index + 1) & 63) : int2(n, n);
    int step = (int)floor(randTime);
    int4 v = indices.xyxy + (int)randSeed + int4(step, step, step + 1, step + 1) * 104729;
    v = (v << 13) ^ v;
    float4 values = (float4)((v * (v * v * a + b) + c) & 0x7fffffff) / 2147483647.0;
    float t = frac(randTime), u = frac(x);
    float2 timed = lerp(values.xy, values.zw, t * t * (3.0 - 2.0 * t));
    return lerp(timed.x, timed.y, u * u * (3.0 - 2.0 * u));
}
float hash1(int n) { return animatedHash(n, 17389, 611953, 1611623773, max(abs(frequency) * 2.0, 1e-6)); }
float hash2(int n) { return animatedHash(n, 27449, 746773, 1824261409, 0.0); }

float randomTurn(float turn, bool inverse) {
    float amount = saturate(randWidth * 0.01);
    if (amount == 0.0) return turn;
    float count = max(abs(frequency), 1.0);
    float shift = (hash2(-19) - 0.5) * amount;
    if (inverse) turn -= shift;
    float x = frac(turn), total = 0.0, accumulated = 0.0;
    [loop] for (int i = 0; i < (int)ceil(count); ++i) {
        float size = min(count - i, 1.0);
        float weight = -log(max(hash2(i + 101), 1e-6)) * size;
        total += weight;
        accumulated += weight * saturate((x * count - i) / size);
    }
    if (!inverse) return floor(turn) + shift + lerp(x, accumulated / total, amount);
    accumulated = 0.0;
    [loop] for (int j = 0; j < (int)ceil(count); ++j) {
        float size = min(count - j, 1.0);
        float weight = -log(max(hash2(j + 101), 1e-6)) * size;
        float width = lerp(size / count, weight / total, amount);
        if (x <= accumulated + width)
            return floor(turn) + (j + size * (x - accumulated) / width) / count;
        accumulated += width;
    }
    return floor(turn) + 1.0;
}

float2 waveLayout(int waveType) {
    float steps = (waveType == WAVE_SAWTOOTH_POS || waveType == WAVE_SAWTOOTH_NEG) ? 1.0 : 2.0;
    float offset = waveType == WAVE_SEMICIRCLE_POS ? 0.5
        : (steps == 1.0 || waveType == WAVE_SEMICIRCLE_NEG) ? 0.0 : 0.25;
    return float2(offset, steps);
}

float2 circularWaveInput(float turns, float2 layout) {
    if (frequency == 0.0) return float2(layout.x, 0.0);
    float count = abs(frequency) * layout.y;
    float turn = randomTurn(turns, true);
    float cycle = floor(turn + 0.5);
    float index = floor((turn - cycle) * count);
    float start = randomTurn(cycle + index / count, false);
    float finish = randomTurn(cycle + (index + 1.0) / count, false);
    float t = saturate((turns - start) / (finish - start));
    return float2(layout.x + (index + t) / layout.y, TAU * (finish - start));
}

float applyWaveHeightRandomness(float wave, float input, float amount, int waveType) {
    if (amount == 0.0) return wave;
    switch (waveType) {
        case WAVE_SQUARE: {
            float amp = 1.0 - hash1((int)floor(2 * input)) * amount;
            return wave * amp;
        }
        case WAVE_SAWTOOTH_POS:
        case WAVE_SAWTOOTH_NEG: {
            float amp = 1.0 - hash1(2 * (int)floor(input) + 1) * amount;
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
            float polarity = waveType == WAVE_SEMICIRCLE_NEG ? -1.0 : 1.0;
            wave *= polarity;
            int shift = waveType == WAVE_SEMICIRCLE_NEG ? 1 : 0;
            float amp0 = 0, amp1 = 0;
            if (frac(input) < 0.5) {
                amp0 = 1.0 - hash1((int)floor(2 * input - 1.0) + shift) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input) + shift) * amount;
                return polarity * ((amp1 + amp0) * wave + (amp1 - amp0)) * 0.5;
            } else {
                amp0 = 1.0 - hash1((int)floor(2 * input - 1.0) + shift) * amount;
                amp1 = 1.0 - hash1((int)floor(2 * input) + shift) * amount;
                return polarity * ((amp1 + amp0) * wave + (amp0 - amp1)) * 0.5;
            }
        }
        default:
            return wave;
    }
}

float polygonRadius(float outer, float inner, float span, float t) {
    if (outer < inner) {
        float swap = outer;
        outer = inner;
        inner = swap;
        t = 1.0 - t;
    }
    if (t <= 0.0) return outer;
    if (t >= 1.0 || outer == inner) return inner;
    float limit = min(span, acos(saturate(inner / outer)));
    if (limit <= 1e-6) return lerp(outer, inner, t);
    float angle = t * span;
    if (angle >= limit) return inner;
    return clamp(inner * outer * sin(limit)
        / (inner * sin(limit - angle) + outer * sin(angle)), inner, outer);
}

float calculatePolygonRadius(float wave, float2 input, int waveType) {
    float inner = max(baseRadius + amplitude * (basePosition
        + applyWaveHeightRandomness(-1.0, input.x, randHeight * 0.01, waveType)), 1e-6);
    float outer = max(baseRadius + amplitude * (basePosition
        + applyWaveHeightRandomness(1.0, input.x, randHeight * 0.01, waveType)), inner);
    return polygonRadius(outer, inner, input.y, saturate((1.0 - wave) * 0.5));
}

float sampleWaveRadius(float turns, int waveType, float2 layout) {
    float2 input = circularWaveInput(turns, layout);
    float wave = waveType == WAVE_SAWTOOTH_POS && frac(input.x) == 0.0
        ? 1.0 : calculateWave(input.x, waveType);
    return calculatePolygonRadius(wave, input, waveType);
}

float circularSawRadius(float turns, int waveType) {
    float radius = sampleWaveRadius(turns, waveType, waveLayout(waveType));
    float count = abs(frequency), whole = floor(count), fraction = frac(count);
    if (fraction == 0.0) return radius;
    float turn = randomTurn(turns, true);
    float cycle = floor(turn + 0.5);
    turn -= cycle;
    float index = floor(whole * 0.5);
    float edge = index / count;
    if (abs(turn) <= edge) return radius;
    if (turn < 0.0) cycle -= 1.0;
    float start = randomTurn(cycle + edge, false);
    float finish = randomTurn(cycle + 1.0 - edge, false);
    float t = saturate((turns - start) / (finish - start));
    float2 input = float2(index + t, TAU * (finish - start));
    if (frac(whole * 0.5) == 0.0)
        return calculatePolygonRadius(calculateWave(input.x, waveType), input, waveType);
    float polarity = waveType == WAVE_SAWTOOTH_POS ? 1.0 : -1.0;
    float2 leftInput = float2(index, 0.0), rightInput = float2(-index - 1.0, 0.0);
    float left = calculatePolygonRadius(-polarity, leftInput, waveType);
    float right = calculatePolygonRadius(polarity, rightInput, waveType);
    float center = randomTurn(cycle + 0.5, false);
    float middle = polygonRadius(left, right, TAU * (finish - start), (center - start) / (finish - start));
    if (turns < center) {
        float tip = lerp(middle, calculatePolygonRadius(polarity, leftInput, waveType), fraction);
        return polygonRadius(left, tip, TAU * (center - start), (turns - start) / (center - start));
    }
    float tip = lerp(middle, calculatePolygonRadius(-polarity, rightInput, waveType), fraction);
    return polygonRadius(tip, right, TAU * (finish - center), (turns - center) / (finish - center));
}

float circularRadius(float turns, int waveType) {
    float2 layout = waveLayout(waveType);
    float radius = sampleWaveRadius(turns, waveType, layout);
    float count = abs(frequency);
    float whole = floor(count), fraction = frac(count);
    if (fraction == 0.0) return radius;
    float turn = randomTurn(turns, true);
    float cycle = floor(turn + 0.5);
    turn -= cycle;
    float edge = whole / (2.0 * count);
    if (abs(turn) <= edge) return radius;
    float side = turn < 0.0 ? -1.0 : 1.0;
    float start = randomTurn(cycle + side * edge, false);
    float finish = randomTurn(cycle + side * 0.5, false);
    float t = saturate((turns - start) / (finish - start));
    float2 startInput = circularWaveInput(start, layout);
    startInput.x = layout.x + side * whole * 0.5;
    float startWave = calculateWave(startInput.x, waveType);
    float startRadius = calculatePolygonRadius(startWave, startInput, waveType);
    float2 endInput = float2(layout.x + count * 0.5, TAU / (count * layout.y));
    if (waveType == WAVE_SQUARE) {
        endInput.x = layout.x + side * (whole + 1.0) * 0.5;
        float endRadius = calculatePolygonRadius(-startWave, endInput, waveType);
        return t < 0.5 ? startRadius : endRadius;
    }
    float endRadius = calculatePolygonRadius(calculateWave(endInput.x, waveType), endInput, waveType);
    float wave = calculateWave(layout.x + (whole + t) * 0.5, waveType);
    float polarity = frac(whole * 0.5) == 0.0 ? 1.0 : -1.0;
    float progress = saturate((1.0 - polarity * wave) * 0.5);
    return polygonRadius(startRadius, endRadius, TAU * abs(finish - start), progress);
}

float4 psmain(float4 pos : SV_Position) : SV_Target {
    float2 p = pos.xy - originalCenter;
    float theta = atan2(p.y, p.x);
    int waveType = (int)(shape + 0.5);
    float turns = frac(theta / TAU + 0.25 + phase + 0.5) - 0.5;
    float radius;
    if (waveType == WAVE_SAWTOOTH_POS || waveType == WAVE_SAWTOOTH_NEG) {
        radius = circularSawRadius(turns, waveType);
    } else {
        radius = circularRadius(turns, waveType);
    }
    float inverseScale = baseRadius / radius;
    float2 sourcePosition = p * inverseScale;
    float2 uv = (originalCenter + sourcePosition) / resolution;
    if (any(uv < 0.0) || any(uv > 1.0)) return float4(0, 0, 0, 0);
    return src.Sample(s, uv);
}
