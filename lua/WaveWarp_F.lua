-- ローカル関数定義
local function createRotationMatrix(angle)
    local cos_a = math.cos(angle)
    local sin_a = math.sin(angle)
    return { cos_a, -sin_a, sin_a, cos_a }
end

local function calculateExpansionParams(height, rot_rad, fix, mirror, center, orig_w, orig_h)
    local cx, cy = tonumber(center[1]) or 0, tonumber(center[2]) or 0
    
    -- 固定処理タイプ別の拡張パラメータ
    local baseX = height * math.abs(math.sin(rot_rad))
    local baseY = height * math.abs(math.cos(rot_rad))
    local base = { top = baseY, bottom = baseY, right = baseX, left = baseX }

    local function map_from_vals(t)
        local map = {
            [1] = {t.top, t.bottom, t.right, t.left}, -- なし(全方向)
            [2] = {0, 0, 0, 0},                       -- すべて
            [3] = {0, 0, t.right, t.left},            -- 上下
            [4] = {t.top, t.bottom, 0, 0},            -- 左右
            [5] = {0, t.bottom, t.right, t.left},     -- 上
            [6] = {t.top, 0, t.right, t.left},        -- 下
            [7] = {t.top, t.bottom, t.right, 0},      -- 左
            [8] = {t.top, t.bottom, 0, t.left},       -- 右
        }
        return map[fix] or {t.top, t.bottom, t.right, t.left}
    end

    if not mirror or mirror <= 0 then
        return map_from_vals(base)
    end

    local r = -rot_rad
    local function proj_abs(x, y)
        return math.abs(x * math.sin(r) + y * math.cos(r))
    end

    -- 共通の分母を計算（最大の投影距離）
    local maxProjection = math.max(
        proj_abs(orig_w + 2 * math.abs(cx), orig_h + 2 * math.abs(cy)),
        proj_abs(orig_w + 2 * math.abs(cx), -orig_h - 2 * math.abs(cy))
    )

    -- 各辺の投影（分子）
    local topNum = math.max(
        proj_abs(orig_w + 2 * cx, orig_h + 2 * cy),
        proj_abs(-orig_w + 2 * cx, orig_h + 2 * cy)
    )
    local bottomNum = math.max(
        proj_abs(orig_w - 2 * cx, orig_h - 2 * cy),
        proj_abs(-orig_w - 2 * cx, orig_h - 2 * cy)
    )
    local rightNum = math.max(
        proj_abs(orig_w - 2 * cx, orig_h - 2 * cy),
        proj_abs(orig_w - 2 * cx, -orig_h - 2 * cy)
    )
    local leftNum = math.max(
        proj_abs(orig_w + 2 * cx, orig_h + 2 * cy),
        proj_abs(orig_w + 2 * cx, -orig_h + 2 * cy)
    )

    -- 比率とスケーリング
    local topRatio    = topNum    / maxProjection
    local bottomRatio = bottomNum / maxProjection
    local rightRatio  = rightNum  / maxProjection
    local leftRatio   = leftNum   / maxProjection

    local scaled = {
        top = base.top * topRatio,
        bottom = base.bottom * bottomRatio,
        right = base.right * rightRatio,
        left = base.left * leftRatio,
    }

    return map_from_vals(scaled)
end

return {
    createRotationMatrix = createRotationMatrix,
    calculateExpansionParams = calculateExpansionParams,
};