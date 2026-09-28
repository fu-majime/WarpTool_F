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

local function encode_values(values)
    local text = {}
    for i = 1, #values do text[i] = string.format("%.6f", values[i]) end
    return table.concat(text, ",")
end

local function decode_values(text)
    local values = {}
    for value in text:gmatch("[^,]+") do values[#values + 1] = tonumber(value) end
    return values
end

local function get_bridge(obj)
    local ok, bridge = pcall(function() return obj.module("ScriptEditBridge") end)
    if ok and type(bridge) == "table" and type(bridge.request) == "function" then
        return bridge
    end
end

-- obj/global は呼び出し元の描画コンテキストを渡す
-- items: {表示名, values内の位置, 最小値, 最大値}
local function syncControls(obj, global, options)
    if not obj.getoption("gui") or options.injected then return end
    local script_name = obj.getoption("script_name")
    local effect_index = 0
    for relative = -1, -64, -1 do
        local ok, name = pcall(function()
            return obj.getoption("script_name", relative, false)
        end)
        if ok and name == script_name then effect_index = effect_index + 1 end
    end
    local key = options.key .. script_name .. ":" .. tostring(effect_index)
    local state = global[key]
    local layout = options.layout
    local from_layout = options.from_layout
    if options.use_anchor then
        if state == "0" or (type(state) == "string" and state:match("^%d,")) then
            local values = state == "0" and options.values or decode_values(state:sub(3))
            state = "l," .. encode_values(options.to_layout(values))
            global[key] = state
        end
        if type(state) == "string" and (state:sub(1, 2) == "l," or state:sub(1, 2) == "w,") then
            local payload = state:sub(3)
            local pending = decode_values(payload)
            local matched = true
            for i = 1, 6 do
                if math.abs(layout[i] - pending[i]) >= 0.00001 then matched = false end
            end
            if not matched then
                if state:sub(1, 2) == "l," then
                    local bridge = get_bridge(obj)
                    if bridge then
                        global[key] = "w," .. payload
                        local called, accepted = pcall(bridge.request, script_name, effect_index, "アンカー制御", payload)
                        if not called or not accepted then global[key] = state end
                    else
                        global[key] = "a," .. encode_values(from_layout(layout))
                        return
                    end
                end
                return nil, pending
            end
        end
        global[key] = "a," .. encode_values(from_layout(layout))
        return
    end
    if type(state) ~= "string" or state == "0" then
        global[key] = "0"
        return
    end
    if state:sub(1, 2) == "l," or state:sub(1, 2) == "w," then
        state = "a," .. encode_values(from_layout(decode_values(state:sub(3))))
    end
    if state:sub(1, 2) == "a," then
        local values = decode_values(state:sub(3))
        for _, item in ipairs(options.items) do
            if item[3] then
                local index = item[2]
                local value = values[index] or options.values[index]
                values[index] = math.max(item[3], math.min(item[4], math.floor(value * 100 + 0.5) / 100))
            end
        end
        state = "1," .. encode_values(values)
        global[key] = state
    end
    local stage, payload = state:match("^(%d),(.+)$")
    stage = tonumber(stage)
    if not stage then return end
    local values = decode_values(payload)
    while stage <= #options.items do
        local item = options.items[stage]
        local index = item[2]
        local changed, value
        if item[3] then
            changed = math.abs(options.values[index] - values[index]) >= 0.005
            value = tostring(values[index])
        else
            changed = math.abs(options.values[index] - values[index]) >= 0.00001
                or math.abs(options.values[index + 1] - values[index + 1]) >= 0.00001
            value = encode_values({values[index], values[index + 1]})
        end
        local next_state = stage == #options.items and "0" or tostring(stage + 1) .. "," .. payload
        if changed then
            local bridge = get_bridge(obj)
            if not bridge then
                global[key] = "0"
                return
            end
            global[key] = next_state
            local called, accepted = pcall(bridge.request, script_name, effect_index, item[1], value)
            if not called or not accepted then global[key] = tostring(stage) .. "," .. payload end
            return values
        end
        global[key] = next_state
        stage = stage + 1
    end
end

local function resolveLayout(layout, injected)
    local result = {}
    injected = type(injected) == "table" and injected or {}
    for i = 1, 6 do result[i] = tonumber(injected[i]) or layout[i] end
    return result
end

return {
    syncControls = syncControls,
    resolveLayout = resolveLayout,
    createRotationMatrix = createRotationMatrix,
    calculateExpansionParams = calculateExpansionParams,
};
