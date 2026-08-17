use std::fmt::Write;

pub const MIN_PINS: usize = 0;
pub const MAX_PINS: usize = 20;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum PinKind {
    #[default]
    Position = 0,
    Bone = 1,
    Bend = 2,
    Detail = 3,
    Starch = 4,
    Overlap = 5,
}

impl PinKind {
    pub const ALL: [Self; 6] = [
        Self::Position,
        Self::Bone,
        Self::Bend,
        Self::Detail,
        Self::Starch,
        Self::Overlap,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::Position => "位置ピン",
            Self::Bone => "ボーンピン",
            Self::Bend => "ベンドピン",
            Self::Detail => "詳細ピン",
            Self::Starch => "スターチピン",
            Self::Overlap => "重なりピン",
        }
    }

    fn from_number(value: usize) -> Option<Self> {
        Self::ALL.into_iter().find(|kind| *kind as usize == value)
    }

    fn has_destination(self) -> bool {
        !matches!(self, Self::Bend | Self::Starch | Self::Overlap)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PinModel {
    pub src: Vec<Point>,
    pub dst: Vec<Point>,
    pub kinds: Vec<PinKind>,
    /// Parent indices for bone pins. Indices are zero-based in Rust and are
    /// converted to the script's one-based nested forest when serialized.
    pub bone_parents: Vec<Option<usize>>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PinOperation {
    Add {
        point: Point,
        kind: PinKind,
        bone_parent: Option<usize>,
    },
    InsertBone {
        point: Point,
        parent: usize,
        child: usize,
    },
    AttachRoot {
        parent: usize,
        root: usize,
    },
    Delete(usize),
    MoveSourceAndDestination {
        index: usize,
        point: Point,
    },
}

impl PinModel {
    pub fn parse(
        pin_count: &str,
        src: &str,
        dst: &str,
        pin_types: &str,
        bone_forest: &str,
    ) -> Result<Self, String> {
        let pin_count = pin_count
            .trim()
            .parse::<usize>()
            .map_err(|_| "pins is not an integer".to_owned())?;
        if !(MIN_PINS..=MAX_PINS).contains(&pin_count) {
            return Err(format!("pins must be between {MIN_PINS} and {MAX_PINS}"));
        }

        let src = parse_points(src, "src")?;
        if src.len() != pin_count {
            return Err(format!(
                "src contains {} pins, expected {pin_count}",
                src.len()
            ));
        }
        let kinds = parse_pin_kinds(pin_types, pin_count);
        let dst = expand_destinations(parse_points(dst, "dst")?, &src, &kinds)?;
        let bone_parents = parse_bone_forest(bone_forest, &kinds);
        Ok(Self {
            src,
            dst,
            kinds,
            bone_parents,
        })
    }

    pub fn parse_or_default(
        pin_count: &str,
        src: &str,
        dst: &str,
        pin_types: &str,
        bone_forest: &str,
        image_size: Option<(usize, usize)>,
    ) -> Result<(Self, bool), String> {
        let pin_count = pin_count
            .trim()
            .parse::<usize>()
            .map_err(|_| "pins is not an integer".to_owned())?;
        if !(MIN_PINS..=MAX_PINS).contains(&pin_count) {
            return Err(format!("pins must be between {MIN_PINS} and {MAX_PINS}"));
        }
        if pin_count == 0 {
            return Self::parse("0", src, dst, pin_types, bone_forest).map(|model| (model, false));
        }
        let kinds = parse_pin_kinds(pin_types, pin_count);
        let destination_count = kinds.iter().filter(|kind| kind.has_destination()).count();
        if !src.trim().is_empty() && (!dst.trim().is_empty() || destination_count == 0) {
            return Self::parse(&pin_count.to_string(), src, dst, pin_types, bone_forest)
                .map(|model| (model, false));
        }
        let (width, height) =
            image_size.ok_or_else(|| "src/dst is empty; waiting for the input image".to_owned())?;

        let source_was_default = src.trim().is_empty();
        let src = if source_was_default {
            default_points(pin_count, width as f32, height as f32)
        } else {
            let src = parse_points(src, "src")?;
            if src.len() != pin_count {
                return Err(format!(
                    "src contains {} pins, expected {pin_count}",
                    src.len()
                ));
            }
            src
        };
        let destination_was_default = dst.trim().is_empty() && destination_count != 0;
        let dst = if destination_was_default {
            src.clone()
        } else {
            expand_destinations(parse_points(dst, "dst")?, &src, &kinds)?
        };
        let bone_parents = parse_bone_forest(bone_forest, &kinds);
        Ok((
            Self {
                src,
                dst,
                kinds,
                bone_parents,
            },
            source_was_default || destination_was_default,
        ))
    }

    pub fn apply(&mut self, operation: PinOperation) -> Result<(), String> {
        match operation {
            PinOperation::Add {
                point,
                kind,
                bone_parent,
            } => {
                validate_point(point)?;
                if self.src.len() >= MAX_PINS {
                    return Err(format!("pin limit ({MAX_PINS}) reached"));
                }
                if kind != PinKind::Bone && bone_parent.is_some() {
                    return Err("only bone pins can have a parent".to_owned());
                }
                if let Some(parent) = bone_parent
                    && self.kinds.get(parent) != Some(&PinKind::Bone)
                {
                    return Err("bone parent is missing or is not a bone pin".to_owned());
                }
                self.src.push(point);
                self.dst.push(point);
                self.kinds.push(kind);
                self.bone_parents.push(bone_parent);
            }
            PinOperation::InsertBone {
                point,
                parent,
                child,
            } => {
                validate_point(point)?;
                if self.src.len() >= MAX_PINS {
                    return Err(format!("pin limit ({MAX_PINS}) reached"));
                }
                if self.kinds.get(parent) != Some(&PinKind::Bone)
                    || self.kinds.get(child) != Some(&PinKind::Bone)
                    || self.bone_parents.get(child) != Some(&Some(parent))
                {
                    return Err("bone connection no longer exists".to_owned());
                }
                let source_parent = self.src[parent];
                let source_child = self.src[child];
                let source_dx = source_child.x - source_parent.x;
                let source_dy = source_child.y - source_parent.y;
                let source_length_squared = source_dx * source_dx + source_dy * source_dy;
                if source_length_squared <= f32::EPSILON {
                    return Err("cannot split a zero-length bone connection".to_owned());
                }
                let projection = (((point.x - source_parent.x) * source_dx
                    + (point.y - source_parent.y) * source_dy)
                    / source_length_squared)
                    .clamp(0.0, 1.0);
                let destination_parent = self.dst[parent];
                let destination_child = self.dst[child];
                let destination = Point {
                    x: destination_parent.x
                        + (destination_child.x - destination_parent.x) * projection,
                    y: destination_parent.y
                        + (destination_child.y - destination_parent.y) * projection,
                };
                validate_point(destination)?;
                let inserted = self.src.len();
                self.src.push(point);
                self.dst.push(destination);
                self.kinds.push(PinKind::Bone);
                self.bone_parents.push(Some(parent));
                self.bone_parents[child] = Some(inserted);
            }
            PinOperation::AttachRoot { parent, root } => {
                if parent == root
                    || self.kinds.get(parent) != Some(&PinKind::Bone)
                    || self.kinds.get(root) != Some(&PinKind::Bone)
                    || self.bone_parents.get(root) != Some(&None)
                {
                    return Err("only another root bone can be attached".to_owned());
                }
                let mut ancestor = Some(parent);
                while let Some(pin) = ancestor {
                    if pin == root {
                        return Err("attaching this root would create a cycle".to_owned());
                    }
                    ancestor = self.bone_parents[pin];
                }
                self.bone_parents[root] = Some(parent);
            }
            PinOperation::Delete(index) => {
                if self.src.len() == MIN_PINS {
                    return Err("there is no pin to delete".to_owned());
                }
                if index >= self.src.len() {
                    return Err("pin index is out of range".to_owned());
                }
                let promoted_parent = self.bone_parents[index]
                    .map(|parent| if parent > index { parent - 1 } else { parent });
                self.src.remove(index);
                self.dst.remove(index);
                self.kinds.remove(index);
                self.bone_parents.remove(index);
                for parent in &mut self.bone_parents {
                    *parent = match *parent {
                        Some(value) if value == index => promoted_parent,
                        Some(value) if value > index => Some(value - 1),
                        value => value,
                    };
                }
            }
            PinOperation::MoveSourceAndDestination { index, point } => {
                validate_point(point)?;
                let Some(source) = self.src.get(index).copied() else {
                    return Err("pin index is out of range".to_owned());
                };
                let Some(destination) = self.dst.get(index).copied() else {
                    return Err("corresponding destination pin is missing".to_owned());
                };
                let moved_destination = Point {
                    x: destination.x + point.x - source.x,
                    y: destination.y + point.y - source.y,
                };
                validate_point(moved_destination)?;
                self.src[index] = point;
                self.dst[index] = moved_destination;
            }
        }
        Ok(())
    }

    pub fn pin_count_string(&self) -> String {
        self.src.len().to_string()
    }

    pub fn src_string(&self) -> String {
        format_points(&self.src)
    }

    pub fn dst_string(&self) -> String {
        format_points(
            &self
                .dst
                .iter()
                .copied()
                .zip(self.kinds.iter().copied())
                .filter_map(|(point, kind)| kind.has_destination().then_some(point))
                .collect::<Vec<_>>(),
        )
    }

    pub fn pin_types_string(&self) -> String {
        self.kinds
            .iter()
            .map(|kind| (*kind as u8).to_string())
            .collect::<Vec<_>>()
            .join(",")
    }

    pub fn bone_forest_string(&self) -> String {
        let mut children = vec![Vec::new(); self.src.len()];
        for (child, parent) in self.bone_parents.iter().copied().enumerate() {
            if let Some(parent) = parent {
                children[parent].push(child);
            }
        }

        fn write_tree(output: &mut String, pin: usize, children: &[Vec<usize>]) {
            output.push('{');
            write!(output, "{}", pin + 1).unwrap();
            for &child in &children[pin] {
                output.push(',');
                write_tree(output, child, children);
            }
            output.push('}');
        }

        let mut output = String::from("{");
        let mut first = true;
        for pin in 0..self.src.len() {
            if self.kinds[pin] == PinKind::Bone && self.bone_parents[pin].is_none() {
                if !first {
                    output.push(',');
                }
                first = false;
                write_tree(&mut output, pin, &children);
            }
        }
        output.push('}');
        output
    }
}

fn parse_pin_kinds(value: &str, count: usize) -> Vec<PinKind> {
    let values = value
        .split([',', '{', '}'])
        .filter_map(|part| part.trim().parse::<usize>().ok())
        .map(|value| PinKind::from_number(value).unwrap_or_default())
        .collect::<Vec<_>>();
    (0..count)
        .map(|index| values.get(index).copied().unwrap_or_default())
        .collect()
}

fn expand_destinations(
    destinations: Vec<Point>,
    sources: &[Point],
    kinds: &[PinKind],
) -> Result<Vec<Point>, String> {
    // Accept the former full-width representation so existing projects migrate
    // losslessly the next time the editor writes the value.
    if destinations.len() == sources.len() {
        return Ok(destinations);
    }

    let expected = kinds.iter().filter(|kind| kind.has_destination()).count();
    if destinations.len() != expected {
        return Err(format!(
            "dst contains {} pins, expected {expected}",
            destinations.len()
        ));
    }

    let mut compact = destinations.into_iter();
    Ok(sources
        .iter()
        .copied()
        .zip(kinds.iter().copied())
        .map(|(source, kind)| {
            if kind.has_destination() {
                compact.next().unwrap()
            } else {
                source
            }
        })
        .collect())
}

#[derive(Debug)]
struct BoneTree {
    pin: usize,
    children: Vec<BoneTree>,
}

fn parse_bone_forest(value: &str, kinds: &[PinKind]) -> Vec<Option<usize>> {
    fn skip_separators(bytes: &[u8], cursor: &mut usize) {
        while bytes
            .get(*cursor)
            .is_some_and(|byte| byte.is_ascii_whitespace() || *byte == b',')
        {
            *cursor += 1;
        }
    }

    fn parse_tree(bytes: &[u8], cursor: &mut usize) -> Option<BoneTree> {
        skip_separators(bytes, cursor);
        if bytes.get(*cursor) != Some(&b'{') {
            return None;
        }
        *cursor += 1;
        skip_separators(bytes, cursor);
        let start = *cursor;
        while bytes.get(*cursor).is_some_and(u8::is_ascii_digit) {
            *cursor += 1;
        }
        let pin = std::str::from_utf8(&bytes[start..*cursor])
            .ok()?
            .parse::<usize>()
            .ok()?
            .checked_sub(1)?;
        let mut children = Vec::new();
        loop {
            skip_separators(bytes, cursor);
            match bytes.get(*cursor) {
                Some(b'{') => children.push(parse_tree(bytes, cursor)?),
                Some(b'}') => {
                    *cursor += 1;
                    return Some(BoneTree { pin, children });
                }
                _ => return None,
            }
        }
    }

    fn visit(
        tree: &BoneTree,
        parent: Option<usize>,
        kinds: &[PinKind],
        parents: &mut [Option<usize>],
        seen: &mut [bool],
    ) {
        if tree.pin >= kinds.len() || kinds[tree.pin] != PinKind::Bone || seen[tree.pin] {
            return;
        }
        seen[tree.pin] = true;
        parents[tree.pin] = parent;
        for child in &tree.children {
            visit(child, Some(tree.pin), kinds, parents, seen);
        }
    }

    let mut parents = vec![None; kinds.len()];
    let mut seen = vec![false; kinds.len()];
    let bytes = value.trim().as_bytes();
    if bytes.first() != Some(&b'{') {
        return parents;
    }
    let mut cursor = 1;
    loop {
        skip_separators(bytes, &mut cursor);
        match bytes.get(cursor) {
            Some(b'{') => {
                let Some(tree) = parse_tree(bytes, &mut cursor) else {
                    return vec![None; kinds.len()];
                };
                visit(&tree, None, kinds, &mut parents, &mut seen);
            }
            Some(b'}') | None => break,
            _ => return vec![None; kinds.len()],
        }
    }
    parents
}

fn default_points(count: usize, width: f32, height: f32) -> Vec<Point> {
    let half_width = width * 0.5;
    let half_height = height * 0.5;
    match count {
        0 => Vec::new(),
        1 => vec![Point::default()],
        2 => vec![
            Point {
                x: -half_width * 0.3,
                y: 0.0,
            },
            Point {
                x: half_width * 0.3,
                y: 0.0,
            },
        ],
        _ => {
            let mut points = Vec::with_capacity(count);
            points.push(Point::default());
            let radius_x = half_width * 0.4;
            let radius_y = half_height * 0.4;
            for index in 0..count - 1 {
                let angle = index as f32 / (count - 1) as f32 * std::f32::consts::TAU
                    - std::f32::consts::FRAC_PI_2;
                points.push(Point {
                    x: angle.cos() * radius_x,
                    y: angle.sin() * radius_y,
                });
            }
            points
        }
    }
}

fn parse_points(value: &str, name: &str) -> Result<Vec<Point>, String> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }
    let values = trimmed
        .split(',')
        .map(|part| {
            part.trim()
                .parse::<f32>()
                .map_err(|_| format!("{name} contains an invalid number"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if values.len() % 2 != 0 {
        return Err(format!("{name} must contain x,y pairs"));
    }

    values
        .chunks_exact(2)
        .map(|pair| {
            let point = Point {
                x: pair[0],
                y: pair[1],
            };
            validate_point(point)?;
            Ok(point)
        })
        .collect()
}

fn validate_point(point: Point) -> Result<(), String> {
    if point.x.is_finite() && point.y.is_finite() {
        Ok(())
    } else {
        Err("pin coordinates must be finite".to_owned())
    }
}

fn format_points(points: &[Point]) -> String {
    let mut output = String::new();
    for (index, point) in points.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        write!(
            &mut output,
            "{},{}",
            format_number(point.x),
            format_number(point.y)
        )
        .unwrap();
    }
    output
}

pub(crate) fn format_number(value: f32) -> String {
    let mut value = format!("{value:.4}");
    while value.contains('.') && value.ends_with('0') {
        value.pop();
    }
    if value.ends_with('.') {
        value.pop();
    }
    if value == "-0" {
        value = "0".to_owned();
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_formats_pin_values() {
        let model = PinModel::parse("2", "-10.5,20,30,-0", "-10.5,20,35,4", "", "").unwrap();
        assert_eq!(model.src[1], Point { x: 30.0, y: -0.0 });
        assert_eq!(model.src_string(), "-10.5,20,30,0");
    }

    #[test]
    fn rejects_invalid_shapes_and_numbers() {
        assert!(PinModel::parse("2", "0,0", "0,0,1,1", "", "").is_err());
        assert!(PinModel::parse("1", "NaN,0", "0,0", "", "").is_err());
        assert!(PinModel::parse("-1", "", "", "", "").is_err());
    }

    #[test]
    fn add_and_delete_keep_source_and_destination_aligned() {
        let mut model = PinModel::parse("1", "0,0", "10,20", "", "").unwrap();
        model
            .apply(PinOperation::Add {
                point: Point { x: 5.0, y: 6.0 },
                kind: PinKind::Detail,
                bone_parent: None,
            })
            .unwrap();
        assert_eq!(model.src[1], model.dst[1]);
        model.apply(PinOperation::Delete(0)).unwrap();
        assert_eq!(model.src, vec![Point { x: 5.0, y: 6.0 }]);
        model.apply(PinOperation::Delete(0)).unwrap();
        assert!(model.src.is_empty());
        assert!(model.apply(PinOperation::Delete(0)).is_err());
    }

    #[test]
    fn zero_pins_parse_without_an_input_image() {
        let (model, synthesized) = PinModel::parse_or_default("0", "", "", "", "", None).unwrap();
        assert!(model.src.is_empty());
        assert!(model.dst.is_empty());
        assert!(model.kinds.is_empty());
        assert!(!synthesized);
    }

    #[test]
    fn moving_source_translates_destination_by_the_same_delta() {
        let mut model = PinModel::parse("1", "0,0", "10,20", "", "").unwrap();
        model
            .apply(PinOperation::MoveSourceAndDestination {
                index: 0,
                point: Point { x: 3.0, y: 4.0 },
            })
            .unwrap();
        assert_eq!(model.src[0], Point { x: 3.0, y: 4.0 });
        assert_eq!(model.dst[0], Point { x: 13.0, y: 24.0 });
    }

    #[test]
    fn enforces_the_twenty_pin_limit() {
        let values = (0..MAX_PINS)
            .flat_map(|index| [index.to_string(), "0".to_owned()])
            .collect::<Vec<_>>()
            .join(",");
        let mut model = PinModel::parse(&MAX_PINS.to_string(), &values, &values, "", "").unwrap();
        assert!(
            model
                .apply(PinOperation::Add {
                    point: Point::default(),
                    kind: PinKind::Position,
                    bone_parent: None,
                })
                .is_err()
        );
    }

    #[test]
    fn creates_the_same_default_layout_as_the_script() {
        let (model, synthesized) =
            PinModel::parse_or_default("3", "", "", "", "", Some((100, 200))).unwrap();
        assert!(synthesized);
        assert_eq!(model.src[0], Point::default());
        assert!((model.src[1].x - 0.0).abs() < 0.001);
        assert!((model.src[1].y + 40.0).abs() < 0.001);
        assert_eq!(model.src, model.dst);
    }

    #[test]
    fn pin_types_and_bone_forest_round_trip() {
        let model = PinModel::parse(
            "6",
            "0,0,1,0,2,0,3,0,4,0,5,0",
            "0,0,1,0,2,0,3,0,4,0,5,0",
            "1,1,2,1,4,5",
            "{{1,{2,{4}}},{3}}",
        )
        .unwrap();
        assert_eq!(
            model.kinds,
            vec![
                PinKind::Bone,
                PinKind::Bone,
                PinKind::Bend,
                PinKind::Bone,
                PinKind::Starch,
                PinKind::Overlap,
            ]
        );
        assert_eq!(
            model.bone_parents,
            vec![None, Some(0), None, Some(1), None, None]
        );
        assert_eq!(model.pin_types_string(), "1,1,2,1,4,5");
        assert_eq!(model.bone_forest_string(), "{{1,{2,{4}}}}");
    }

    #[test]
    fn destination_serialization_skips_non_movable_kinds_without_losing_ids() {
        let model = PinModel::parse(
            "6",
            "0,0,1,0,0,1,-1,0,0,-1,-1,1",
            "2,3,4,-2,-2,-3",
            "0,3,4,0,5,2",
            "",
        )
        .unwrap();

        assert_eq!(model.dst.len(), 6);
        assert_eq!(model.dst[0], Point { x: 2.0, y: 3.0 });
        assert_eq!(model.dst[1], Point { x: 4.0, y: -2.0 });
        assert_eq!(model.dst[2], model.src[2]);
        assert_eq!(model.dst[3], Point { x: -2.0, y: -3.0 });
        assert_eq!(model.dst[4], model.src[4]);
        assert_eq!(model.dst[5], model.src[5]);
        assert_eq!(model.dst_string(), "2,3,4,-2,-2,-3");
    }

    #[test]
    fn full_width_legacy_destinations_are_migrated_when_serialized() {
        let model = PinModel::parse("3", "0,0,1,0,2,0", "10,0,11,0,12,0", "0,2,3", "").unwrap();

        assert_eq!(model.dst_string(), "10,0,12,0");
    }

    #[test]
    fn deleting_a_bone_repairs_parent_ids_and_promotes_its_children() {
        let mut model = PinModel::parse(
            "4",
            "0,0,1,0,2,0,3,0",
            "0,0,1,0,2,0,3,0",
            "1,1,1,0",
            "{{1,{2,{3}}}}",
        )
        .unwrap();
        model.apply(PinOperation::Delete(1)).unwrap();
        assert_eq!(
            model.kinds,
            vec![PinKind::Bone, PinKind::Bone, PinKind::Position]
        );
        assert_eq!(model.bone_parents, vec![None, Some(0), None]);
        assert_eq!(model.bone_forest_string(), "{{1,{2}}}");
    }

    #[test]
    fn inserting_a_bone_splits_the_selected_connection() {
        let mut model =
            PinModel::parse("2", "0,0,100,0", "10,20,210,20", "1,1", "{{1,{2}}}").unwrap();
        model
            .apply(PinOperation::InsertBone {
                point: Point { x: 40.0, y: 0.0 },
                parent: 0,
                child: 1,
            })
            .unwrap();
        assert_eq!(model.kinds, vec![PinKind::Bone; 3]);
        assert_eq!(model.bone_parents, vec![None, Some(2), Some(0)]);
        assert_eq!(model.dst[2], Point { x: 90.0, y: 20.0 });
        assert_eq!(model.bone_forest_string(), "{{1,{3,{2}}}}");
    }

    #[test]
    fn attaching_a_root_moves_its_entire_subtree_and_rejects_cycles() {
        let mut model = PinModel::parse(
            "5",
            "0,0,1,0,2,0,3,0,4,0",
            "0,0,1,0,2,0,3,0,4,0",
            "1,1,1,1,1",
            "{{1,{2}},{3,{4,{5}}}}",
        )
        .unwrap();
        model
            .apply(PinOperation::AttachRoot { parent: 1, root: 2 })
            .unwrap();
        assert_eq!(model.bone_forest_string(), "{{1,{2,{3,{4,{5}}}}}}");
        assert!(
            model
                .apply(PinOperation::AttachRoot { parent: 4, root: 0 })
                .is_err()
        );
        assert!(
            model
                .apply(PinOperation::AttachRoot { parent: 0, root: 3 })
                .is_err()
        );
    }
}
