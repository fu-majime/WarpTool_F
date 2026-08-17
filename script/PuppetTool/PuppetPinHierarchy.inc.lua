-- Hierarchical pin solver.
--
-- This file does not access obj or global.  The including script owns the UI,
-- persistence, and application of the solved positions.
local pin_hierarchy = {}

local EPSILON = 0.1

function pin_hierarchy.default_positions(count)
    local out = {}
    for pin = 1, count do
        out[pin * 2 - 1] = (pin - (count + 1) * 0.5) * 80
        out[pin * 2] = 0
    end
    return out
end

local function point(values, pin)
    local i = pin * 2 - 1
    return tonumber(values[i]), tonumber(values[i + 1])
end

local function distance(ax, ay, bx, by)
    local dx, dy = ax - bx, ay - by
    return math.sqrt(dx * dx + dy * dy)
end

local function atan2(y, x)
    if x > 0 then return math.atan(y / x) end
    if x < 0 then return math.atan(y / x) + (y >= 0 and math.pi or -math.pi) end
    if y > 0 then return math.pi * 0.5 end
    if y < 0 then return -math.pi * 0.5 end
    return 0
end

local function normalize(angle)
    while angle > math.pi do angle = angle - math.pi * 2 end
    while angle <= -math.pi do angle = angle + math.pi * 2 end
    return angle
end

-- Convert a nested forest into arrays. Invalid, duplicate, and cyclic
-- entries are ignored. Pins omitted from the forest become additional roots.
local function compile_forest(value, count)
    local parent, children, order, seen = {}, {}, {}, {}
    for pin = 1, count do children[pin] = {} end

    local function walk(tree, parent_pin)
        if type(tree) ~= "table" then return end
        local pin = math.floor(tonumber(tree[1]) or 0)
        if pin < 1 or pin > count or seen[pin] then return end
        seen[pin], parent[pin] = true, parent_pin
        order[#order + 1] = pin
        if parent_pin then
            children[parent_pin][#children[parent_pin] + 1] = pin
        end
        for i = 2, #tree do walk(tree[i], pin) end
    end

    if type(value) == "table" then
        for i = 1, #value do walk(value[i], nil) end
    end
    for pin = 1, count do
        if not seen[pin] then
            seen[pin] = true
            order[#order + 1] = pin
        end
    end
    return parent, children, order
end

local function hierarchy_signature(parent, count)
    local signature = count * 1009
    for pin = 1, count do signature = signature + pin * ((parent[pin] or 0) + 31) end
    return signature
end

-- State contains one local angle and one wait flag per pin. A signature
-- invalidates it when the pin count or forest changes.
local function load_state(serialized, count, signature, parent, order, positions)
    local values = {}
    if type(serialized) == "string" then
        for token in string.gmatch(serialized, "[^,]+") do
            values[#values + 1] = tonumber(token)
        end
    end
    if #values == count * 2 + 3 and values[1] == 1 and
       values[2] == count and values[3] == signature then
        local angles, waiting = {}, {}
        for pin = 1, count do
            angles[pin] = values[3 + pin]
            waiting[pin] = values[3 + count + pin] ~= 0
        end
        return angles, waiting
    end

    local angles, waiting, global_angle = {}, {}, {}
    for _, pin in ipairs(order) do
        local parent_pin = parent[pin]
        if parent_pin then
            local x, y = point(positions, pin)
            local px, py = point(positions, parent_pin)
            local absolute = atan2(y - py, x - px)
            angles[pin] = normalize(absolute - (global_angle[parent_pin] or 0))
            global_angle[pin] = absolute
        else
            angles[pin], global_angle[pin] = 0, 0
        end
        waiting[pin] = false
    end
    return angles, waiting
end

local function save_state(count, signature, angles, waiting)
    local values = {1, count, signature}
    for pin = 1, count do values[#values + 1] = angles[pin] or 0 end
    for pin = 1, count do values[#values + 1] = waiting[pin] and 1 or 0 end
    local out = {}
    for i = 1, #values do out[i] = string.format("%.17g", values[i]) end
    return table.concat(out, ",")
end

function pin_hierarchy.solve(source, destination, hierarchy, count, serialized, editing_source)
    local parent, children, order = compile_forest(hierarchy, count)
    local signature = hierarchy_signature(parent, count)
    local angles, waiting = load_state(serialized, count, signature, parent, order, destination)
    local lengths = {}
    for pin = 1, count do
        local parent_pin = parent[pin]
        if parent_pin then
            local x, y = point(source, pin)
            local px, py = point(source, parent_pin)
            lengths[pin] = distance(x, y, px, py)
        end
    end

    local function pose()
        local x, y, global_angle = {}, {}, {}
        for _, pin in ipairs(order) do
            local parent_pin = parent[pin]
            if parent_pin then
                global_angle[pin] = (global_angle[parent_pin] or 0) + angles[pin]
                x[pin] = x[parent_pin] + lengths[pin] * math.cos(global_angle[pin])
                y[pin] = y[parent_pin] + lengths[pin] * math.sin(global_angle[pin])
            else
                x[pin], y[pin] = point(source, pin)
                global_angle[pin] = 0
            end
        end
        return x, y, global_angle
    end

    local x, y, global_angle = pose()

    local function subtree_arrived(pin)
        local dx, dy = point(destination, pin)
        if distance(dx, dy, x[pin], y[pin]) >= EPSILON then return false end
        for _, child in ipairs(children[pin]) do
            if not subtree_arrived(child) then return false end
        end
        return true
    end

    for pin = 1, count do
        if waiting[pin] and subtree_arrived(pin) then waiting[pin] = false end
    end

    if not editing_source then
        -- Parents occur before children. The first changed node wins, and a
        -- waiting ancestor shields its descendants from stale async values.
        for _, pin in ipairs(order) do
            local parent_pin = parent[pin]
            if parent_pin then
                local blocked = false
                local ancestor = parent_pin
                while ancestor do
                    if waiting[ancestor] then blocked = true; break end
                    ancestor = parent[ancestor]
                end
                local dx, dy = point(destination, pin)
                if not blocked and distance(dx, dy, x[pin], y[pin]) >= EPSILON then
                    local absolute = atan2(dy - y[parent_pin], dx - x[parent_pin])
                    angles[pin] = normalize(absolute - global_angle[parent_pin])
                    waiting[pin] = #children[pin] > 0
                    break
                end
            end
        end
        x, y = pose()
    end

    local solved = {}
    for pin = 1, count do
        solved[pin * 2 - 1], solved[pin * 2] = x[pin], y[pin]
    end
    return solved, save_state(count, signature, angles, waiting)
end

function pin_hierarchy.encode(values)
    local out = {}
    for i = 1, #values do
        local value = math.abs(values[i]) < EPSILON * 0.5 and 0 or values[i]
        out[i] = string.format("%.4f", value)
    end
    return table.concat(out, ",")
end

function pin_hierarchy.changed(current, solved)
    for i = 1, #solved do
        if not tonumber(current[i]) or math.abs(current[i] - solved[i]) >= EPSILON then
            return true
        end
    end
    return false
end

