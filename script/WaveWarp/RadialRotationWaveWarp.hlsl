Texture2D src : register(t0);
cbuffer constant0 : register(b0) {
    float2 resolution;
    float2 originalResolution;
    float2 center;
    float2 originalCenter;
    float time;
    float period;
    float radAmp;
    float radFreq;
    float shape;
    float radPhase;
    float spiralMix;
    float rotAmp;
    float rotRadius;
    float rotPhase;
    float debug;
    float debug_A;
    float debug_B;
};

SamplerState s;

static const float TAU = 6.28318530718;

static const int WAVE_SINE = 1;
static const int WAVE_SQUARE = 2;
static const int WAVE_TRIANGLE = 3;
static const int WAVE_SAWTOOTH_POS = 4;
static const int WAVE_SAWTOOTH_NEG = 5;
static const int WAVE_CIRCLE = 6;
static const int WAVE_SEMICIRCLE_POS = 7;
static const int WAVE_SEMICIRCLE_NEG = 8;

float inverse_x_plus_log(float y);
float inverse_weighted_log_linear(float y, float alpha, float beta);

float hash1(int n) {
    n += 12345; // simple seed
    n = (n << 13) ^ n;
    return (float)((n * (n * n * 15731 + 789221) + 1376312627) & 0x7fffffff) / 2147483647.0;
}

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
    float2 texCoord = pos.xy / resolution;
    
    float2 centerUV = originalCenter / resolution;
    float2 p = (texCoord - centerUV) * resolution;
    
    float r = length(p);
    float theta = atan2(p.y, p.x);

    
    float waveInput = (theta / TAU + radPhase) * radFreq;
    if (abs(period) > 0) {
        waveInput += time / period * radFreq;
    }
    
    int waveType = (int)(shape + 0.5);
    float wave = calculateWave(waveInput, waveType);
    
    float rotInput = r / rotRadius + rotPhase;
    if (abs(period) > 0) rotInput += time / period;
    float rotWave = calculateWave(rotInput, waveType);

    float2 newTexCoord = texCoord;
    // 放射モード
    if (false)
    {
        float n = max(1.0, abs(radFreq));
        float BaseRadius = max(originalResolution.x/2, originalResolution.y/2);
        float Rmax = max(BaseRadius + radAmp, 1e-6);
        float Rmin = max(BaseRadius - radAmp, 1e-6);
        float pi_over_n = TAU / 2 / n;
        float tan_alpha = (Rmax / Rmin - cos(pi_over_n)) / sin(pi_over_n);
        float alpha = atan(tan_alpha);
        float term = TAU * (wave + 1.0) / 4.0;
        
        float limit_ratio = (1.0 - cos(pi_over_n)) / (1.0 + cos(pi_over_n));
        float h_limit = BaseRadius * limit_ratio;
        float straightness = smoothstep(0.0, h_limit, abs(radAmp));
        straightness = clamp(abs(radAmp) / h_limit, 0.0, 1.0);
        // straightness = 1;

        float denom_poly = cos(term / n - alpha);
        float denom_circle = cos(alpha);
        float denom = lerp(denom_circle, denom_poly, straightness);
        float newR = Rmax * cos(alpha) / denom;
        float2 newP = p * (BaseRadius / (BaseRadius + radAmp) / 2.0) / newR;
        float2 aspectCorrection = float2((2.0 * BaseRadius + 2.0 * radAmp) / resolution.x, (2.0 * BaseRadius + 2.0 * radAmp) / resolution.y);
        newTexCoord = centerUV + newP * aspectCorrection;
    
        if (debug > 0.5) {
            float inner = max(BaseRadius - abs(radAmp), 1e-6);
            float outer = max(BaseRadius + abs(radAmp), inner + 1e-6);
            float ginfo = smoothstep(inner, outer, r); // 0=中心側、1=外側（radAmp影響域の割合）
            return float4(frac(radFreq * theta / TAU), 0, -1 * 0.5 + 0.5, 1.0);
        }
    }
    // 回転モード
    else if (false)
    {
        float angleOffset = rotWave * rotAmp * TAU;
        float newTheta = theta - angleOffset;
        float2 newP_rot = float2(cos(newTheta), sin(newTheta)) * r;
        newTexCoord = centerUV + newP_rot / resolution;
    }
    // 螺旋モード
    {
        // スパイラル座標系に変換
        float twistSin = sin(spiralMix * TAU / 4.0);
        float twistCos = cos(spiralMix * TAU / 4.0);
        float twistTan = tan(spiralMix * TAU / 4.0);

        // 極座標の計算（rはlogをとる）
        float logR = log(r + 1e-6);

        // 基本のu, vを計算
        float u = twistCos * logR - twistSin * (theta + 0*TAU * round((logR / twistTan - theta) / TAU));
        float v = twistSin * logR + twistCos * (theta + 0*TAU * round((logR / twistTan - theta) / TAU));
        

        // 波形を計算
        float phase_adjusted = radPhase;
        float rot_phase_adjusted = rotPhase;
        if (abs(period) > 0) {
            phase_adjusted += time / period;
            rot_phase_adjusted -= time / period;
        }
        float waveInput = inverse_weighted_log_linear(
            v + TAU * twistCos * (phase_adjusted - rot_phase_adjusted / radFreq),
            1 - twistCos, TAU * twistCos / radFreq / rotRadius) / rotRadius + rot_phase_adjusted;
        if (abs(twistCos) > 0.01) {
            waveInput = frac(radFreq * v / twistCos / TAU);
        } else {
            waveInput = exp(v) / rotRadius + rotPhase;
        }
        float rotWaveInput = (twistSin * r + twistCos) / rotRadius + rot_phase_adjusted;
        // waveInput = (v / TAU + radPhase) * radFreq;
        // waveInput = exp(v) / rotRadius + rotPhase;

        float wave = calculateWave(waveInput, waveType);
        float rotWave = calculateWave(rotWaveInput, waveType);

        float n = max(1.0, abs(radFreq));
        float BaseRadius = max(originalResolution.x/2, originalResolution.y/2);
        float Rmax = max(BaseRadius + radAmp, 1e-6);
        float Rmin = max(BaseRadius - radAmp, 1e-6);
        float pi_over_n = TAU / 2 / n;
        float tan_alpha = (Rmax / Rmin - cos(pi_over_n)) / sin(pi_over_n);
        float alpha = atan(tan_alpha);
        float term = TAU * (wave + 1.0) / 4.0;
        
        float limit_ratio = (1.0 - cos(pi_over_n)) / (1.0 + cos(pi_over_n));
        float h_limit = BaseRadius * limit_ratio;
        float straightness = smoothstep(0.0, h_limit, abs(radAmp));
        straightness = clamp(abs(radAmp) / h_limit, 0.0, 1.0);
        // straightness = 1;

        float denom_poly = cos(term / n - alpha);
        float denom_circle = cos(alpha);
        float denom = lerp(denom_circle, denom_poly, straightness);
        float newR = Rmax * cos(alpha) / denom;
        float angleOffset = wave * rotAmp * TAU;
        
        float newU = u;
        newU += (1 - twistSin) * log(BaseRadius/ (Rmax * cos(alpha)) * denom );
        newU += twistSin * (wave * rotAmp * TAU) + 0 * rotWave * rotAmp * TAU;
        // newU += calculateWave(radFreq * v / twistCos / TAU, waveType);
        float2 newP = float2(
            exp(twistCos * newU + twistSin * v) * cos(-twistSin * newU +twistCos * v),
            exp(twistCos * newU + twistSin * v) * sin(-twistSin * newU +twistCos * v)
        );
        newTexCoord = centerUV + newP / resolution;

        if (debug > 0.5) {
            float u_comp = (abs(twistSin) > 0.01) ? (u / twistSin) : u;
            float newu_comp = (abs(twistSin) > 0.01) ? (newU / twistSin) : newU;
            float v_comp = (abs(twistCos) > 0.01) ? (v / twistCos) : v;

            return float4(
                1 * frac(v_comp / TAU),
                0.2 * frac(radFreq * newu_comp / TAU),
                frac(radFreq * u_comp / TAU),
                1.0
            ) * 0.8 + src.Sample(s, newTexCoord) * 0.2;
        }
    }
    
    if (debug > 0.5) {
        return float4(frac(theta / TAU), 0.0, wave * 0.5 + 0.5, 1.0);
    }

    if (newTexCoord.x < 0.0 || newTexCoord.x > 1.0 || 
        newTexCoord.y < 0.0 || newTexCoord.y > 1.0) 
    {
        return float4(0, 0, 0, 0); // 範囲外は透明
    }
    return src.Sample(s, newTexCoord);
}

// x + ln(x) = y を満たす正の x を求める
float inverse_x_plus_log(float y) {
    // 初期値: 大きい y では漸近展開、小さめなら単純な近似を使う
    float x = 0.0;
    if (y > 8.0) {
        float logy = log(y);
        x = y - logy + logy / max(y, 1e-3); // 1次漸近 + 補正
    } else if (y > 1.0) {
        x = y - log(y);
    } else {
        x = exp(y); // y が小さい（負を含む）領域では x ≈ e^y
    }

    // 数値安定のために下限を入れる
    x = max(x, 1e-6);

    // ニュートン法で 3-4 回程度改良（収束が速いので固定回数）
    [unroll]
    for (int i = 0; i < 4; ++i) {
        float f = x + log(x) - y;
        float fp = 1.0 + 1.0 / x;
        float dx = f / fp;
        x -= dx;
        x = max(x, 1e-6);
    }

    return x;
}
// alpha*ln(x) + beta*x = y を満たす正の x を求める
float inverse_weighted_log_linear(float y, float alpha, float beta) {
    if (abs(beta) < 1e-6) {
        return exp(y / alpha);
    }
    if (abs(alpha) < 1e-6) {
        return y / beta;
    }

    float V = (y / alpha) - (log(alpha) - log(beta));
    float w = inverse_x_plus_log(V);
    return (alpha / beta) * w;
}
