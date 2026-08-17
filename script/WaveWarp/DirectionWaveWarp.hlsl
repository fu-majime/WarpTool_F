Texture2D src : register(t0);
cbuffer constant0 : register(b0) {
    float2 resolution;
    float2 originalResolution;
    float2 center;
    float2 originalCenter;
    float time;
    float height;
    float width;
    float rot;
    float2 rot1_row0;
    float2 rot1_row1;
    float2 rot2_row0;
    float2 rot2_row1;
    float period;
    float shape;
    float phase;
    float mirror;
    float fix;
    float fadePixel;
    float randHeight;
    float randWidth;
    float randSeed;
    float debug;
};

SamplerState s;

// 定数定義
static const float TAU = 6.28318530718;
static const float TOLLERANCE = 1e-6;

// 固定タイプ定数
static const int FIX_NONE = 1;
static const int FIX_EDGE_ALL = 2;
static const int FIX_EDGE_VERTICAL = 3;
static const int FIX_EDGE_HORIZONTAL = 4;
static const int FIX_EDGE_TOP = 5;
static const int FIX_EDGE_BOTTOM = 6;
static const int FIX_EDGE_LEFT = 7;
static const int FIX_EDGE_RIGHT = 8;
static const int FIX_CENTER = 9;

// 波形タイプ定数
static const int WAVE_SINE = 1;
static const int WAVE_SQUARE = 2;
static const int WAVE_TRIANGLE = 3;
static const int WAVE_SAWTOOTH_POS = 4;
static const int WAVE_SAWTOOTH_NEG = 5;
static const int WAVE_CIRCLE = 6;
static const int WAVE_SEMICIRCLE_POS = 7;
static const int WAVE_SEMICIRCLE_NEG = 8;

// 波のモード定数
static const int WAVE_MODE_NORMAL = 0;
static const int WAVE_MODE_MIRROR = 1;

// 関数宣言
float linearstep(float edge0, float edge1, float x);
float bezierTransition(float x, float slope, float bend);
float smoothTransition(float x, float slope);
float emphasizedTransition(float x, float slope);

float2 applyRotation(float2 coord, float2x2 rotMatrix);

float2 calculateFade(float2 texCoord, float2 warpCoord, float fadePixel, int edgeType);
float2 applyEdgeClamp(float2 coord, int edgeType);

float calculateWaveInput(float2 rotatedCoord, float width, float phase, float period, float time);
float applyWaveWidthRandomness(float input, float amount, int waveType);
float calculateWave(float input, int waveType);
float applyWaveHeightRandomness(float wave, float input, float amount, int waveType);
float2 calculateWarpedCoord(float2 texCoord, float2 rotatedCoord,
                            float wave, float height, float width, float rot, float2x2 rot2,
                            int waveMode, int fixType);
                            
// 整数ハッシュ関数 (0.0 ~ 1.0 を返す)
float hash1(int n) {
    n += (int)randSeed;
    n = (n << 13) ^ n;
    return (float)((n * (n * n * 17389 + 611953) + 1611623773) & 0x7fffffff) / 2147483647.0;
}
float hash2(int n) {
    n += (int)randSeed;
    n = (n << 13) ^ n;
    return (float)((n * (n * n * 27449 + 746773) + 1824261409) & 0x7fffffff) / 2147483647.0;
}

float4 psmain(float4 pos : SV_Position) : SV_Target {
    float2 texCoord = pos.xy / resolution;
    
    // 回転行列を構築
    float2x2 rot1 = float2x2(rot1_row0, rot1_row1);
    float2x2 rot2 = float2x2(rot2_row0, rot2_row1);
    
    // 1. 最初の回転を適用
    float2 rotatedCoord = applyRotation(texCoord, rot1);
    
    // 2. 波形入力値を計算
    float waveInput = calculateWaveInput(rotatedCoord, width, phase, period, time);
    waveInput = applyWaveWidthRandomness(waveInput, randWidth * 0.01, (int)(shape + 0.5));
    
    // 3. 波形を計算
    int waveType = (int)(shape + 0.5);
    float wave = calculateWave(waveInput, waveType);
    wave = applyWaveHeightRandomness(wave, waveInput, randHeight * 0.01, waveType);
    
    // 4. ワープされた座標を計算
    int edgeType = (int)(fix + 0.5);
    int waveMode = (mirror > 0.5) ? WAVE_MODE_MIRROR : WAVE_MODE_NORMAL;
    float2 finalCoord = calculateWarpedCoord(texCoord, rotatedCoord, wave, height, width, rot, rot2, waveMode, edgeType);
    
    // デバッグモード分岐
    if (debug > 0.5) {
    }
    
    return src.Sample(s, finalCoord);
}

// 回転を適用する関数
float2 applyRotation(float2 coord, float2x2 rotMatrix) {
    // originalCenterを正規化座標系での回転中心として使用
    float2 centerNormalized = originalCenter / resolution;
    return mul((coord - centerNormalized) * resolution, rotMatrix) / resolution + centerNormalized;
}

///// calculateFade用ヘルパー関数
// 1.0 - bezierTransition(1.0 - x, 1.0, 1 - slope/3)
// 数学的には: x + x² - x³
float simpleTransitionSlope1Inverted(float x) {
    x = clamp(x, 0.0, 1.0);
    float x2 = x * x;
    float x3 = x2 * x;
    return x + x2 - x3;
}

float reinhard_weight(float x, float knee_threshold) {
    return (knee_threshold + max(0.0, -x)) / (knee_threshold + max(0.0, x));
}

// エッジフェード係数を計算するヘルパー関数
float2 computeEdgeFade(float normalizedDist, float simpleDist, float warpAmount, float fadeWidth) {
    // float fadeEdge = 1.0 - bezierTransition(1.0 - simpleDist, 1.0, 2.0/3.0);
    float fadeEdge = simpleTransitionSlope1Inverted(simpleDist);
    float slope = reinhard_weight(warpAmount, fadeWidth);
    float fadeWarp = emphasizedTransition(normalizedDist, slope);
    return float2(fadeWarp, fadeEdge);
}
///// calculateFade用ヘルパー関数

float2 calculateFade(float2 texCoord, float2 warpCoord, float fadePixel, int edgeType) {
    float fadeWidthU = fadePixel / resolution.x;
    float fadeWidthV = fadePixel / resolution.y;

    // 左エッジ: X軸の正方向からのフェード (warp, edge)
    float leftNormalizedDist = texCoord.x / (fadeWidthU + max(0.0, -warpCoord.x));
    float leftSimpleDist = texCoord.x / fadeWidthU;
    float2 fadeLeft = computeEdgeFade(leftNormalizedDist, leftSimpleDist, warpCoord.x, fadeWidthU);

    // 右エッジ: X軸の負方向からのフェード (warp, edge)
    float rightNormalizedDist = (1.0 - texCoord.x) / (fadeWidthU + max(0.0, warpCoord.x));
    float rightSimpleDist = (1.0 - texCoord.x) / fadeWidthU;
    float2 fadeRight = computeEdgeFade(rightNormalizedDist, rightSimpleDist, -warpCoord.x, fadeWidthU);

    // 上エッジ: Y軸の正方向からのフェード (warp, edge)
    float topNormalizedDist = texCoord.y / (fadeWidthV + max(0.0, -warpCoord.y));
    float topSimpleDist = texCoord.y / fadeWidthV;
    float2 fadeTopRaw = computeEdgeFade(topNormalizedDist, topSimpleDist, warpCoord.y, fadeWidthV);
    float2 fadeTop = float2(fadeTopRaw.y, fadeTopRaw.x);

    // 下エッジ: Y軸の負方向からのフェード (warp, edge)
    float bottomNormalizedDist = (1.0 - texCoord.y) / (fadeWidthV + max(0.0, warpCoord.y));
    float bottomSimpleDist = (1.0 - texCoord.y) / fadeWidthV;
    float2 fadeBottomRaw = computeEdgeFade(bottomNormalizedDist, bottomSimpleDist, -warpCoord.y, fadeWidthV);
    float2 fadeBottom = float2(fadeBottomRaw.y, fadeBottomRaw.x);

    // 中央固定: 中央のとき0.0、一番端で1.0にする
    // 中央の定義は、端のfadePixel分を除いた範囲
    float2 fadeCenter;
    {
        float centerMinX = fadeWidthU + warpCoord.x;
        float centerMaxX = 1.0 - fadeWidthU + warpCoord.x;
        float centerMinY = fadeWidthV + warpCoord.y;
        float centerMaxY = 1.0 - fadeWidthV + warpCoord.y;
        if (texCoord.x >= centerMinX && texCoord.x <= centerMaxX &&
            texCoord.y >= centerMinY && texCoord.y <= centerMaxY) {
            fadeCenter = float2(0.0, 0.0);
        } else {
            fadeCenter = float2(1.0, 1.0);
        }
    }
    {
        float centerMinX = fadeWidthU;
        float centerMaxX = 1.0 - fadeWidthU;
        float centerMinY = fadeWidthV;
        float centerMaxY = 1.0 - fadeWidthV;
        
        // X軸方向のフェード計算
        float fadeCenterX;
        if (texCoord.x >= centerMinX && texCoord.x <= centerMaxX) {
            fadeCenterX = 0.0; // 中央領域内は0.0
        } else if (texCoord.x < centerMinX) {
            // 左端からのフェード (端で1.0、中央境界で0.0)
            fadeCenterX = 1.0 - simpleTransitionSlope1Inverted(texCoord.x / fadeWidthU);
        } else {
            // 右端からのフェード (端で1.0、中央境界で0.0)
            fadeCenterX = 1.0 - simpleTransitionSlope1Inverted((1.0 - texCoord.x) / fadeWidthU);
        }
        
        // Y軸方向のフェード計算
        float fadeCenterY;
        if (texCoord.y >= centerMinY && texCoord.y <= centerMaxY) {
            fadeCenterY = 0.0; // 中央領域内は0.0
        } else if (texCoord.y < centerMinY) {
            // 上端からのフェード (端で1.0、中央境界で0.0)
            fadeCenterY = 1.0 - simpleTransitionSlope1Inverted(texCoord.y / fadeWidthV);
        } else {
            // 下端からのフェード (端で1.0、中央境界で0.0)
            fadeCenterY = 1.0 - simpleTransitionSlope1Inverted((1.0 - texCoord.y) / fadeWidthV);
        }
        
        // 両軸の最大値を取る（矩形フェード領域）
        float fadeCenterValue = max(fadeCenterX, fadeCenterY);
        fadeCenter = float2(fadeCenterValue, fadeCenterValue);
    }


  

    switch (edgeType) {
        case FIX_EDGE_ALL:        return min(fadeLeft * fadeRight, fadeTop * fadeBottom);
        case FIX_EDGE_VERTICAL:   return fadeTop * fadeBottom;
        case FIX_EDGE_HORIZONTAL: return fadeLeft * fadeRight;
        case FIX_EDGE_TOP:        return fadeTop;
        case FIX_EDGE_BOTTOM:     return fadeBottom;
        case FIX_EDGE_LEFT:       return fadeLeft;
        case FIX_EDGE_RIGHT:      return fadeRight;
        case FIX_CENTER:          return fadeCenter;
        default:                  return float2(1.0, 1.0);
    }
}

// エッジタイプに応じた座標クランプ関数
float2 applyEdgeClamp(float2 coord, int edgeType) {
    float minBoundX = 1.0 / resolution.x;  // 最小境界（1ピクセル）
    float maxBoundX = 1.0 - minBoundX;
    float minBoundY = 1.0 / resolution.y;
    float maxBoundY = 1.0 - minBoundY;
    
    switch (edgeType) {
        case FIX_EDGE_ALL:
            // 全方向クランプ
            return float2(clamp(coord.x, minBoundX, maxBoundX), clamp(coord.y, minBoundY, maxBoundY));
            
        case FIX_EDGE_VERTICAL:
            // Y方向のみクランプ
            return float2(coord.x, clamp(coord.y, minBoundY, maxBoundY));
            
        case FIX_EDGE_HORIZONTAL:
            // X方向のみクランプ
            return float2(clamp(coord.x, minBoundX, maxBoundX), coord.y);
            
        case FIX_EDGE_TOP:
            // 上エッジのみクランプ
            return float2(coord.x, max(coord.y, minBoundY));
            
        case FIX_EDGE_BOTTOM:
            // 下エッジのみクランプ
            return float2(coord.x, min(coord.y, maxBoundY));
            
        case FIX_EDGE_LEFT:
            // 左エッジのみクランプ
            return float2(max(coord.x, minBoundX), coord.y);
            
        case FIX_EDGE_RIGHT:
            // 右エッジのみクランプ
            return float2(min(coord.x, maxBoundX), coord.y);
            
        case FIX_NONE:
        default:
            // クランプなし
            return coord;
    }
}

float calculateWaveInput(float2 rotatedCoord, float width, float phase, float period, float time) {
    float baseX = (rotatedCoord.x - originalCenter.x / resolution.x) * resolution.x;
    
    float waveInput = baseX / width + phase;
    if (abs(period) > 0) {
        waveInput += time / period;
    }
    
    return waveInput;
}

float ReLU(float x) {
    return max(0.0, x);
}

// ワープオフセットのみを計算する関数
float calculateWarpOffset(float2 rotatedCoord, float wave, float height, float rot, int waveMode) {
    float baseWarp = height / resolution.y * wave;
    switch (waveMode) {
        case WAVE_MODE_MIRROR: {
            float sinr = abs(sin(rot));
            float cosr = abs(cos(rot));
            // 鏡面反射ワープ処理
            float2 normalizedCenter = originalCenter / resolution;
            float relativeY = 2.0 * (normalizedCenter.y - rotatedCoord.y);

            float rotatedExtentUV = ((originalResolution.x + 2 * abs(center.x)) * sinr + (originalResolution.y + 2 * abs(center.y)) * cosr) / resolution.y;
            float warpOffset = relativeY * baseWarp / (2.0 * baseWarp + rotatedExtentUV);

            return warpOffset;
        }

        case WAVE_MODE_NORMAL:
        default: {
            // 通常のワープ処理
            return baseWarp;
        }
    }
}

// ワープされた座標を計算する関数
float2 calculateWarpedCoord(float2 texCoord, float2 rotatedCoord,
                            float wave, float height, float width, float rot, float2x2 rot2,
                            int waveMode, int fixType) {
    float warpOffset = calculateWarpOffset(rotatedCoord, wave, height, rot, waveMode);
    {
        float A = fadePixel + height * 2;
        // warpOffset = warpOffset * ReLU((A / resolution.y  - texCoord.y) / (warpOffset + (A - height) / resolution.y));
    }
    float2 warpCoord = mul(float2(0, warpOffset) * resolution, rot2) / resolution;
    float2 fade = calculateFade(texCoord, warpCoord, fadePixel, fixType);
    float2 warped = texCoord + warpCoord * fade;
    return applyEdgeClamp(warped, fixType);
}

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
            return sin(TAU * input); // デフォルトはサイン
    }
}

// linearstep関数の実装
// edge0=0, edge1=0+0のとき、[0,1]で1を返す
float linearstep(float edge0, float edge1, float x) {
    if (edge0 == edge1) {
        return x >= edge0 ? 1.0 : 0.0;
    }
    return clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
}

float bezierTransition(float x, float slope, float bend) {
    x = clamp(x, 0.0, 1.0);
    
    // ベジエ曲線の陽関数解: 傾き0からslopeに変化
    // 3*bend*t^2 + (1-3*bend)*t^3 の逆関数
    
    // 判別式の計算
    float discriminant_base = (3.0 - 2.0*slope - 3.0*bend);
    float discriminant = discriminant_base * discriminant_base - 4.0*slope*(-3.0 + slope + 3.0*bend)*x;
    
    // 平方根の計算（負の場合は0にクランプ）
    float sqrt_discriminant = sqrt(max(0.0, discriminant));
    
    // t の計算
    float denominator = 2.0 * (-3.0 + slope + 3.0*bend);
    float t = (-discriminant_base - sqrt_discriminant) / denominator;
    
    // 結果の計算: 3*bend*t^2 + (1-3*bend)*t^3
    return t*t * (3.0*bend + (1.0 - 3.0*bend) * t);
}

// 滑らかな遷移曲線関数
// x=0で値0、x=1で値1,傾きslopeになる滑らかな曲線
// bezierTransition bend = 1/3
float smoothTransition(float x, float slope) {
    x = clamp(x, 0.0, 1.0);  // 境界値処理を統合
    
    if (slope <= 1.0) {
        // 開始時傾き無限大: (2s-1)x + (2-2s)√x
        return (2.0 * slope - 1.0) * x + (2.0 - 2.0 * slope) * sqrt(x);
    }
    
    if (abs(slope - 2.0) < TOLLERANCE) {
        // 特別ケース: x²
        return x * x;
    }
    
    // 開始時傾き0: ((s-1) - √(1 - s(s-2)(x-1))) / (s-2))²
    float t = ((slope - 1.0) - sqrt(max(0.0, 1.0 - slope * (slope - 2.0) * (x - 1.0)))) / (slope - 2.0);
    return t * t;
}

float emphasizedTransition(float x, float slope) {
    // xの範囲をあらかじめ安全にしておく
    x = clamp(x, 0.0, 1.0);
    
    // slope <= 1.0 のときは数学的に一定値 (2√x - x)
    if (slope <= 1.0) {
        return 2.0 * sqrt(x) - x;
    }

    // 通常ケース
    // g(x, a) = ( f(x, a) - a * f(x, 1) ) / (1 - a)
    // ここで f(x, 1) = x であるため、単純に slope * x となる
    return (smoothTransition(x, slope) - slope * x) / (1.0 - slope);
}
