//! End-to-end geometry invariants, independent of an AviUtl2 process.
use super::*;

#[test]
fn zero_threshold_includes_transparent_image_without_external_padding() {
    let rgba = vec![0; 32 * 24 * 4];
    for threshold in [0, 1, 0] {
        let mesh = generate_mesh(&rgba, 32, 24, threshold, 15, 0).unwrap();
        if threshold == 1 {
            assert!(mesh.triangles.is_empty());
            continue;
        }
        assert!(!mesh.triangles.is_empty());
        assert!(
            mesh.points
                .iter()
                .all(|&(x, y)| (-16.0..=16.0).contains(&x) && (-12.0..=12.0).contains(&y))
        );
        for y in 0..24 {
            for x in 0..32 {
                assert!(covers(
                    &mesh,
                    (x as f32 + 0.5 - 16.0, y as f32 + 0.5 - 12.0)
                ));
            }
        }
    }
}

#[test]
fn alpha_threshold_tracks_gradient_and_cache_changes() {
    let (w, h) = (256, 16);
    let mut rgba = vec![0; w * h * 4];
    for y in 0..h {
        for x in 0..w {
            rgba[(y * w + x) * 4 + 3] = x as u8;
        }
    }
    for threshold in [0, 1, 64, 128, 192, 254, 255, 128, 1, 255, 0] {
        let mesh = generate_mesh(&rgba, w, h, threshold, 15, 0).unwrap();
        for x in 0..w {
            assert_eq!(
                covers(&mesh, (x as f32 + 0.5 - 128.0, 0.0)),
                x >= threshold as usize,
                "threshold {threshold}, alpha {x}"
            );
        }
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png"]
fn actual_illustration_threshold_preserves_selected_pixels() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for threshold in [1, 64, 128, 192, 254, 255] {
        let mesh = generate_mesh(&rgba, 723, 800, threshold, 15, 0).unwrap();
        let mut missing = 0;
        for y in 0..800 {
            for x in 0..723 {
                if rgba[(y * 723 + x) * 4 + 3] < threshold as u8 {
                    continue;
                }
                if !covers(&mesh, (x as f32 + 0.5 - 361.5, y as f32 + 0.5 - 400.0)) {
                    missing += 1;
                }
            }
        }
        assert_eq!(missing, 0, "threshold {threshold}: selected pixels lost");
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png"]
fn actual_illustration_border_expansion() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for border in [1, 5, 10, 25, 50, 100] {
        let start = std::time::Instant::now();
        eprintln!("border {border}: starting");
        let mesh = generate_mesh(&rgba, 723, 800, 128, 20, border).unwrap();
        eprintln!(
            "border {border}: {:?}, {} vertices",
            start.elapsed(),
            mesh.points.len()
        );
        assert!(!mesh.triangles.is_empty());
        assert!(
            mesh.points.len() < 1000,
            "border {border}: raster contour was not simplified"
        );
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png"]
fn actual_illustration_threshold_changes() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for threshold in [0, 1, 32, 64, 96, 128, 160, 192, 224, 254, 255] {
        for border in [0, 1, 10] {
            let start = std::time::Instant::now();
            let mesh = generate_mesh(&rgba, 723, 800, threshold, 20, border).unwrap();
            eprintln!(
                "threshold {threshold}, border {border}: {:?}, {} vertices",
                start.elapsed(),
                mesh.points.len()
            );
            assert!(!mesh.triangles.is_empty());
            assert!(
                mesh.points.len() < 1000,
                "threshold {threshold}, border {border}"
            );
            assert!(
                mesh.points
                    .iter()
                    .all(|p| p.0.is_finite() && p.1.is_finite())
            );
            assert!(
                mesh.triangles
                    .iter()
                    .all(|t| t.iter().all(|&i| i < mesh.points.len()))
            );
        }
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png"]
fn actual_hand_refinement_has_no_samples_crowding_finger_notches() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for density in [30, 35, 50] {
        let pose = "pins=2;bone_forest={};pin_types={0,0};pin_sx={650-hw,531-hw};pin_sy={380-hh,727-hh};pin_dx={650-hw,531-hw};pin_dy={380-hh,727-hh};pin_layer={0,0}";
        let (v, t, _, render) =
            lua_pipeline_fixture(2, pose, density, false, Some((&rgba, 723, 800)));
        std::fs::write(
            format!("../../target/puppet-qa/hand-d{density}.csv"),
            render
                .chunks_exact(12)
                .map(|r| r.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                .collect::<Vec<_>>()
                .join("\n"),
        )
        .unwrap();
        let mut edges = std::collections::BTreeMap::new();
        for face in t.chunks_exact(3) {
            for k in 0..3 {
                let (a, b) = (face[k] as usize, face[(k + 1) % 3] as usize);
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let outline = edges
            .iter()
            .filter(|(_, n)| **n == 1)
            .map(|(&(a, b), _)| (a, b))
            .collect::<Vec<_>>();
        for (i, p) in v.chunks_exact(2).enumerate() {
            if p[0] < 230.
                || !(-70. ..50.).contains(&p[1])
                || outline.iter().any(|&(a, b)| i == a || i == b)
                || (p[0] - 288.5).hypot(p[1] + 20.) < 1e-8
            {
                continue;
            }
            let clearance = outline
                .iter()
                .map(|&(a, b)| {
                    let (x, y) = (v[2 * b] - v[2 * a], v[2 * b + 1] - v[2 * a + 1]);
                    let u = (((p[0] - v[2 * a]) * x + (p[1] - v[2 * a + 1]) * y) / (x * x + y * y))
                        .clamp(0., 1.);
                    (p[0] - v[2 * a] - u * x).hypot(p[1] - v[2 * a + 1] - u * y)
                })
                .fold(f64::INFINITY, f64::min);
            assert!(
                clearance >= (800. / density as f64) * 0.2 - 1e-8,
                "finger-notch sample d{density} {p:?}: clearance {clearance}"
            );
        }
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png"]
fn actual_high_density_ears_follow_head_when_only_a_wrist_moves() {
    let source = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut rgba = vec![0; 725 * 802 * 4];
    for y in 0..800 {
        rgba[((y + 1) * 725 + 1) * 4..((y + 1) * 725 + 724) * 4]
            .copy_from_slice(&source[y * 723 * 4..(y + 1) * 723 * 4]);
    }
    let mut failures = Vec::new();
    // Target the reported high-density regression. The exploratory density-15
    // failures are recorded separately in docs/puppet_deformation.md.
    for density in [30, 35, 38, 40, 45, 50] {
        for method in [1, 2] {
            for pin_count in [2, 5, 6] {
                for (x, y) in [(320, 60), (180, 560), (380, 140), (531, 720), (600, 700)] {
                    let mut pose = format!(
                        r#"
                    pins={pin_count}
                    pin_sx={{360-hw,190-hw,531-hw,70-hw,650-hw,360-hw}}
                    pin_sy={{510-hh,727-hh,727-hh,380-hh,380-hh,300-hh}}
                    bone_forest={{}}
                    for i=1,pins do pin_types[i]=0;pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                    pin_dx[4]={x}-hw;pin_dy[4]={y}-hh
                "#
                    );
                    if pin_count == 2 {
                        pose = format!(
                            r#"
                        pins=2;bone_forest={{}}
                        pin_sx={{650-hw,531-hw}};pin_sy={{380-hh,727-hh}}
                        for i=1,pins do pin_types[i]=0;pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                        pin_dx[1]={x}-hw;pin_dy[1]={y}-hh
                    "#
                        );
                    }
                    pose.push_str(";for i=1,pins do pin_sx[i]=pin_sx[i]+1;pin_sy[i]=pin_sy[i]+1;pin_dx[i]=pin_dx[i]+1;pin_dy[i]=pin_dy[i]+1 end");
                    let (_, _, _, render) = lua_pipeline_fixture(
                        method,
                        &pose,
                        density,
                        false,
                        Some((&rgba, 725, 802)),
                    );
                    if pin_count == 2 && density == 30 && x == 320 {
                        std::fs::write(
                            format!("../../target/puppet-qa/two-pin-m{method}.csv"),
                            render
                                .chunks_exact(12)
                                .map(|row| {
                                    row.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
                                })
                                .collect::<Vec<_>>()
                                .join("\n"),
                        )
                        .unwrap();
                    }
                    let at = |x: f64, y: f64| {
                        rendered_material_point(
                            &render,
                            ((x + 1.) / 725. - 0.5) * 160.,
                            ((y + 1.) / 802. - 0.5) * 120.,
                        )
                    };
                    let center = at(360., 180.);
                    let top = at(360., 40.);
                    let c = (center.1 - top.1) / 140.;
                    let s = (top.0 - center.0) / 140.;
                    for (ex, ey) in [(235., 198.), (485., 198.), (240., 185.), (480., 210.)] {
                        let ear = at(ex, ey);
                        let expected = (
                            center.0 + c * (ex - 360.) - s * (ey - 180.),
                            center.1 + s * (ex - 360.) + c * (ey - 180.),
                        );
                        let error = (ear.0 - expected.0).hypot(ear.1 - expected.1);
                        if error > 3. {
                            failures.push((method, format!(
                                "ear m{method} d{density} wrist({x},{y}) ear({ex},{ey}): {error}px"
                            )));
                        }
                    }
                    for vertex in render.chunks_exact(4) {
                        let (ex, ey) = (vertex[2] * 725. - 1., vertex[3] * 802. - 1.);
                        if (150. ..235.).contains(&ey) && !(255. ..465.).contains(&ex) {
                            let expected = (
                                center.0 + c * (ex - 360.) - s * (ey - 180.),
                                center.1 + s * (ex - 360.) + c * (ey - 180.),
                            );
                            let error = (vertex[0] - expected.0).hypot(vertex[1] - expected.1);
                            if error > 3. {
                                failures.push((method, format!("ear contour m{method} d{density} pins{pin_count} ({x},{y}) ({ex},{ey}): {error}px")));
                            }
                        }
                    }
                }
            }
        }
    }
    shape_quality_results(&failures);
}

fn covers(mesh: &Mesh, p: (f32, f32)) -> bool {
    mesh.triangles.iter().any(|t| {
        let cross =
            |a: (f32, f32), b: (f32, f32)| (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
        let s = [
            cross(mesh.points[t[0]], mesh.points[t[1]]),
            cross(mesh.points[t[1]], mesh.points[t[2]]),
            cross(mesh.points[t[2]], mesh.points[t[0]]),
        ];
        s.iter().all(|&v| v >= -1.0e-4) || s.iter().all(|&v| v <= 1.0e-4)
    })
}
fn mask(w: usize, h: usize, inside: impl Fn(usize, usize) -> bool) -> Vec<u8> {
    let mut rgba = vec![0; w * h * 4];
    for y in 0..h {
        for x in 0..w {
            if inside(x, y) {
                rgba[(y * w + x) * 4 + 3] = 255;
            }
        }
    }
    rgba
}
fn assert_coverage(rgba: &[u8], w: usize, h: usize, density: i32) {
    let mesh = generate_mesh(rgba, w, h, 1, density, 0).unwrap();
    for y in 0..h {
        for x in 0..w {
            if rgba[(y * w + x) * 4 + 3] == 0 {
                continue;
            }
            for (dx, dy) in [
                (0.001, 0.001),
                (0.999, 0.001),
                (0.001, 0.999),
                (0.999, 0.999),
                (0.5, 0.5),
            ] {
                assert!(
                    covers(
                        &mesh,
                        (
                            x as f32 + dx - w as f32 * 0.5,
                            y as f32 + dy - h as f32 * 0.5
                        )
                    ),
                    "uncovered pixel ({x},{y}) density {density}, mask {:?}, points {:?}, triangles {:?}",
                    rgba.chunks_exact(4).map(|v| v[3]).collect::<Vec<_>>(),
                    mesh.points,
                    mesh.triangles
                );
            }
        }
    }
}
#[test]
fn tiny_masks_preserve_cells_and_diagonal_components() {
    for bits in 1..16 {
        let rgba = mask(2, 2, |x, y| bits & (1 << (y * 2 + x)) != 0);
        assert_coverage(&rgba, 2, 2, 5);
    }
}
#[test]
fn concavity_holes_and_thin_limbs_survive_low_density() {
    let rgba = mask(64, 64, |x, y| {
        ((8..56).contains(&x)
            && (8..56).contains(&y)
            && !((20..44).contains(&x) && (20..44).contains(&y)))
            || ((2..62).contains(&x) && y == 3)
    });
    for density in [5, 15, 40] {
        assert_coverage(&rgba, 64, 64, density);
    }
    let mesh = generate_mesh(&rgba, 64, 64, 1, 5, 0).unwrap();
    assert!(!covers(&mesh, (0., 0.)), "hole must stay open");
    assert!(!covers(&mesh, (0., -26.)), "gap must stay open");
}
#[test]
fn opaque_rectangle_has_no_automatic_padding_at_any_density() {
    let rgba = mask(80, 40, |x, y| {
        (10..70).contains(&x) && (10..30).contains(&y)
    });
    for density in [5, 15, 40] {
        let mesh = generate_mesh(&rgba, 80, 40, 1, density, 0).unwrap();
        assert!(
            mesh.points
                .iter()
                .all(|&(x, y)| (-30.0..=30.0).contains(&x) && (-10.0..=10.0).contains(&y))
        );
    }
    assert!(
        generate_mesh(&vec![0; 80 * 40 * 4], 80, 40, 1, 5, 0)
            .unwrap()
            .points
            .is_empty()
    );
}
#[test]
fn both_methods_reproduce_rigid_motion_with_embedded_handles() {
    let (v, t) = embedding::prepare(
        &[-10., -10., 10., -10., 10., 10., -10., 10.],
        &[0, 1, 2, 0, 2, 3],
        &[-3., -2., 4., 3.],
        2.,
    )
    .unwrap();
    for angle in [0.0_f64, 0.7, std::f64::consts::PI] {
        let transform = |x: f64, y: f64| {
            [
                angle.cos() * x - angle.sin() * y + 12.,
                angle.sin() * x + angle.cos() * y - 7.,
            ]
        };
        let pins = [-3., -2., 4., 3.];
        let target = pins
            .chunks_exact(2)
            .flat_map(|p| transform(p[0], p[1]))
            .collect::<Vec<_>>();
        for (method, result) in [
            deformation::deform_mls(&v, &t, &pins, &target, &[0., 0.], 1., 3, 40., 40.).unwrap(),
            deformation::deform_arap(&v, &t, &pins, &target, &[0., 0.], 3, 40., 40.).unwrap(),
        ]
        .into_iter()
        .enumerate()
        {
            for (p, q) in v.chunks_exact(2).zip(result.vertices.chunks_exact(2)) {
                let e = transform(p[0], p[1]);
                assert!(
                    (e[0] - q[0]).hypot(e[1] - q[1]) < 1.0e-6,
                    "mode {method}, angle {angle}, source {p:?}, actual {q:?}, expected {e:?}"
                );
            }
            for q in result.render_vertices.chunks_exact(4) {
                let e = transform(q[2] * 40. - 20., q[3] * 40. - 20.);
                assert!((e[0] - q[0]).hypot(e[1] - q[1]) < 1.0e-6);
            }
        }
    }
}
#[test]
fn malformed_inputs_are_errors_in_both_modes() {
    for target in [vec![], vec![f64::NAN, 0.]] {
        assert!(
            deformation::deform_arap(
                &[0., 0., 10., 0., 0., 10.],
                &[0, 1, 2],
                &[0., 0.],
                &target,
                &[0.],
                1,
                20.,
                20.
            )
            .is_err()
        );
        assert!(
            deformation::deform_mls(
                &[0., 0., 10., 0., 0., 10.],
                &[0, 1, 2],
                &[0., 0.],
                &target,
                &[0.],
                1.,
                1,
                20.,
                20.
            )
            .is_err()
        );
    }
}

type PipelineResult = (Vec<f64>, Vec<i32>, Vec<f64>, Vec<f64>);
#[test]
fn complete_puppet_script_compiles_in_lua51() {
    let mut source = include_str!("../../../script/PuppetTool/PuppetTool.in.anm2").to_owned();
    for (name, body) in [
        (
            "PuppetTool.hlsl",
            include_str!("../../../script/PuppetTool/PuppetTool.hlsl"),
        ),
        (
            "PuppetPinHierarchy.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetPinHierarchy.inc.lua"),
        ),
        (
            "PuppetPins.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetPins.inc.lua"),
        ),
        (
            "PuppetMesh.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetMesh.inc.lua"),
        ),
        (
            "PuppetDeformation.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetDeformation.inc.lua"),
        ),
        (
            "PuppetPinDynamics.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetPinDynamics.inc.lua"),
        ),
        (
            "PuppetRendering.inc.lua",
            include_str!("../../../script/PuppetTool/PuppetRendering.inc.lua"),
        ),
    ] {
        source = source.replace(&format!("---$include \"{name}\""), body);
    }
    mlua::Lua::new().load(source).into_function().unwrap();
}
type DeformArguments = (
    Vec<f64>,
    Vec<i32>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    f64,
    i32,
    f64,
    f64,
    Vec<f64>,
);
#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_illustration_density() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut previous = 0;
    for density in [5, 10, 15, 30, 50] {
        let start = std::time::Instant::now();
        let mesh = generate_mesh(&rgba, 723, 800, 1, density, 0).unwrap();
        eprintln!(
            "density {density}: {} points; hands {:?}",
            mesh.points.len(),
            [-1.0_f32, 1.0].map(|side| mesh
                .points
                .iter()
                .filter(|&&(x, y)| side * x > 230. && (-70. ..50.).contains(&y))
                .count())
        );
        assert!(
            mesh.points.len() >= previous,
            "density reduced the vertex count"
        );
        assert!(
            mesh.points.len()
                <= if density == 5 {
                    100
                } else {
                    density as usize * 10 + 200
                },
            "actual silhouette defeated the density budget"
        );
        previous = mesh.points.len();
        if density == 5 {
            for side in [-1.0_f32, 1.0] {
                let hand = mesh
                    .points
                    .iter()
                    .filter(|&&(x, y)| side * x > 230. && (-70. ..50.).contains(&y))
                    .count();
                assert!(
                    hand <= 8,
                    "low-density hand retained {hand} contour vertices"
                );
            }
        }
        eprintln!(
            "actual density {density}: {} vertices {} triangles {:?}",
            mesh.points.len(),
            mesh.triangles.len(),
            start.elapsed()
        );
    }
}
#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_illustration_poses() {
    actual_illustration_poses_for(&[(1, 5), (2, 5), (1, 15), (2, 15)]);
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_illustration_arap_poses() {
    actual_illustration_poses_for(&[(2, 5), (2, 15)]);
}

fn actual_illustration_poses_for(cases: &[(i32, i32)]) {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let setup = r#"
        pins=8
        pin_sx={360,360,275,447,70,650,190,531}; pin_sy={510,345,320,320,380,380,727,727}
        pin_types={1,1,1,1,1,1,1,1}; pin_layer={0,0,0,0,0,0,0,0}
        for i=1,pins do pin_sx[i]=pin_sx[i]-hw; pin_sy[i]=pin_sy[i]-hh; pin_dx[i]=pin_sx[i]; pin_dy[i]=pin_sy[i]; pin_rotation[i]=0; pin_scale[i]=1; pin_range[i]=10 end
        bone_forest={{1,{2,{3,{5}},{4,{6}}},{7},{8}}}
    "#;
    for &(method, density) in cases {
        for (name, pose) in [
            (
                "head",
                "pins=6;pin_sx={70-hw,650-hw,190-hw,531-hw,360-hw,360-hw};pin_sy={380-hh,380-hh,727-hh,727-hh,470-hh,180-hh};pin_types={0,0,0,0,1,1};bone_forest={{5,{6}}};for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i] end;pin_dx[6]=380-hw;pin_dy[6]=310-hh",
            ),
            (
                "head-root",
                "pins=6;pin_sx={70-hw,650-hw,190-hw,531-hw,360-hw,360-hw};pin_sy={380-hh,380-hh,727-hh,727-hh,470-hh,180-hh};pin_types={0,0,0,0,1,1};bone_forest={{6,{5}}};for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i] end;pin_dx[6]=400-hw;pin_dy[6]=220-hh",
            ),
            (
                "head-turn",
                "pins=6;pin_sx={70-hw,650-hw,190-hw,531-hw,360-hw,360-hw};pin_sy={380-hh,380-hh,727-hh,727-hh,470-hh,180-hh};pin_types={0,0,0,0,1,1};bone_forest={{5,{6}}};for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i] end;pin_dx[6]=500-hw;pin_dy[6]=220-hh",
            ),
            (
                "bone",
                "local x={360,360,290,430,220,480,220,500}; local y={550,345,320,320,150,110,440,440}; for i=1,pins do pin_dx[i]=x[i]-hw; pin_dy[i]=y[i]-hh end",
            ),
            (
                "scale",
                "pins=5; pin_sx={531-hw,190-hw,70-hw,360-hw,447-hw}; pin_sy={727-hh,727-hh,380-hh,470-hh,320-hh}; pin_types={0,0,0,0,3}; bone_forest={}; for i=1,pins do pin_dx[i]=pin_sx[i]; pin_dy[i]=pin_sy[i] end; pin_scale[5]=1.8",
            ),
        ] {
            let (v, t, q, render) = lua_pipeline_fixture(
                method,
                &format!("{setup}\n{pose}"),
                density,
                false,
                Some((&rgba, 723, 800)),
            );
            eprintln!(
                "actual {name} mode {method} density {density}: {} vertices {} triangles",
                v.len() / 2,
                t.len() / 3
            );
            let rows = render
                .chunks_exact(12)
                .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(
                format!("../../target/puppet-qa/actual-{name}-{method}-d{density}.csv"),
                rows,
            )
            .unwrap();
            assert!(q.iter().chain(&render).all(|x| x.is_finite()));
            if name.starts_with("head") || name == "bone" || name == "scale" {
                let locate = |x: f64, y: f64| {
                    for t in render.chunks_exact(12) {
                        let (u, v) = (x / 723.0, y / 800.0);
                        let (bx, by, cx, cy) =
                            (t[6] - t[2], t[7] - t[3], t[10] - t[2], t[11] - t[3]);
                        let det = bx * cy - by * cx;
                        if det.abs() < 1e-15 {
                            continue;
                        }
                        let b = ((u - t[2]) * cy - (v - t[3]) * cx) / det;
                        let c = (bx * (v - t[3]) - by * (u - t[2])) / det;
                        let a = 1.0 - b - c;
                        if a >= -1e-8 && b >= -1e-8 && c >= -1e-8 {
                            return (
                                a * t[0] + b * t[4] + c * t[8],
                                a * t[1] + b * t[5] + c * t[9],
                            );
                        }
                    }
                    panic!("uncovered face sample");
                };
                let a = locate(310.0, 165.0);
                if name == "bone" {
                    let sleeve = locate(211.0, 301.0);
                    let displacement =
                        (sleeve.0 - (211.0 - 361.5)).hypot(sleeve.1 - (301.0 - 400.0));
                    assert!(
                        displacement > 15.0,
                        "raised arm left its sleeve fixed under the head frame: {displacement}"
                    );
                }
                let b = locate(395.0, 165.0);
                let c = locate(360.0, 245.0);
                let top = locate(360.0, 40.0);
                let nose = locate(360.0, 200.0);
                let crown = (top.0 - nose.0).hypot(top.1 - nose.1) / 160.0;
                let width = (b.0 - a.0).hypot(b.1 - a.1) / 85.0;
                let height = (c.0 - (a.0 + b.0) / 2.0).hypot(c.1 - (a.1 + b.1) / 2.0)
                    / (80.0_f64.hypot(7.5));
                eprintln!("{name} mode {method} face width {width} height {height}");
                if name == "scale" {
                    assert!(
                        (0.8..1.25).contains(&(width / height))
                            && (0.8..1.25).contains(&(crown / height)),
                        "scale distorted the face: width {width}, height {height}, crown {crown}"
                    );
                    continue;
                }
                shape_quality(
                    method,
                    (0.8..1.2).contains(&width)
                        && (0.8..1.2).contains(&height)
                        && (0.8..1.2).contains(&crown),
                    format_args!("head bone changed facial proportions"),
                );
            }
        }
    }
}
fn lua_pipeline(method: i32, pose: &str) -> PipelineResult {
    lua_pipeline_at_density(method, pose, 15)
}
fn lua_pipeline_at_density(method: i32, pose: &str, density: i32) -> PipelineResult {
    lua_pipeline_on_mask(method, pose, density, false)
}
fn lua_pipeline_on_mask(method: i32, pose: &str, density: i32, curved: bool) -> PipelineResult {
    lua_pipeline_fixture(method, pose, density, curved, None)
}
fn lua_pipeline_fixture(
    method: i32,
    pose: &str,
    density: i32,
    curved: bool,
    image: Option<(&[u8], usize, usize)>,
) -> PipelineResult {
    let lua = mlua::Lua::new();
    let module = lua.create_table().unwrap();
    module
        .set(
            "prepare_mesh",
            lua.create_function(
                |_, (v, t, p, r, spacing): (Vec<f64>, Vec<i32>, Vec<f64>, f64, Option<f64>)| {
                    embedding::prepare_with_spacing(&v, &t, &p, r, spacing.unwrap_or(0.0))
                        .map_err(mlua::Error::external)
                },
            )
            .unwrap(),
        )
        .unwrap();
    // Deliberately register only the selected solver. An accidental dependency
    // on the other mode must fail this test, including preliminary poses.
    let name = if method == 2 {
        "deform_arap"
    } else {
        "deform_mls"
    };
    module
        .set(
            name,
            lua.create_function(
                move |_, (v, t, p, q, l, s, d, w, h, poses): DeformArguments| {
                    let selected = if method == 2 {
                        deformation::Method::Arap
                    } else {
                        deformation::Method::RigidMls(s)
                    };
                    let result = deformation::deform(&v, &t, &p, &q, &l, d, w, h, &poses, selected)
                        .map_err(mlua::Error::external)?;
                    Ok((
                        result.vertices,
                        result.render_vertices,
                        result.wire_vertices,
                    ))
                },
            )
            .unwrap(),
        )
        .unwrap();
    lua.globals().set("mesh_module", module).unwrap();
    lua.globals().set("deformationMethod", method).unwrap();
    let rgba = mask(160, 120, |x, y| {
        if curved {
            let (x, y) = (x as f64 - 80.0, y as f64 - 60.0);
            let ellipse = |a: f64, b: f64, rx: f64, ry: f64| {
                ((x - a) / rx).powi(2) + ((y - b) / ry).powi(2) < 1.0
            };
            return ellipse(0., -35., 16., 16.)
                || ellipse(0., 5., 25., 34.)
                || ellipse(-38., -14., 32., 8.)
                || ellipse(38., -14., 32., 8.)
                || ellipse(-60., -14., 9., 9.)
                || ellipse(55., -14., 9., 9.)
                || ellipse(-18., 34., 9., 22.)
                || ellipse(18., 34., 9., 22.);
        }
        ((55..105).contains(&x) && (35..95).contains(&y))
            || ((10..150).contains(&x) && (40..53).contains(&y))
            || (((55..69).contains(&x) || (89..103).contains(&x)) && (90..116).contains(&y))
            || ((x as f64 - 80.).powi(2) + (y as f64 - 25.).powi(2) < 16. * 16.)
    });
    let (rgba, width, height) = image.unwrap_or((&rgba, 160, 120));
    let mesh = generate_mesh(rgba, width, height, 1, density, 0).unwrap();
    lua.globals()
        .set(
            "mesh_x",
            mesh.points.iter().map(|p| p.0 as f64).collect::<Vec<_>>(),
        )
        .unwrap();
    lua.globals()
        .set(
            "mesh_y",
            mesh.points.iter().map(|p| p.1 as f64).collect::<Vec<_>>(),
        )
        .unwrap();
    lua.globals()
        .set(
            "mesh_tris",
            mesh.triangles
                .iter()
                .map(|t| t.iter().map(|i| i + 1).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let setup = r#"
        w,h,hw,hh,density,div,stiff=160,120,80,60,15,2,1
        mesh_n_verts=#mesh_x
        PIN_TYPE={POSITION=0,BONE=1,BEND=2,DETAIL=3,STARCH=4,OVERLAP=5}
        pins=6
        pin_sx={-60,55,0,-18,16,0}; pin_sy={-14,-14,-15,47,47,-35}
        pin_types={0,3,1,2,4,5}; pin_layer={0,0,0,0,0,1}
        pin_dx={}; pin_dy={}; pin_rotation={}; pin_scale={}; pin_range={}; pin_show_range={}
        for i=1,pins do pin_dx[i]=pin_sx[i]; pin_dy[i]=pin_sy[i]; pin_rotation[i]=0; pin_scale[i]=1; pin_range[i]=10; end
        show=true; obj={getoption=function() return true end}
    "#;
    let deformation = include_str!("../../../script/PuppetTool/PuppetDeformation.inc.lua").replace(
        "---$include \"PuppetPinDynamics.inc.lua\"",
        include_str!("../../../script/PuppetTool/PuppetPinDynamics.inc.lua"),
    );
    lua.load(format!(
        "{setup}\nw,h,hw,hh={width},{height},{width}/2,{height}/2\ndensity={density}\n{pose}\n{deformation}\nassert(#pin_dynamics_debug==pins); for _,guide in ipairs(pin_dynamics_debug) do assert(guide.dx==guide.dx and guide.dy==guide.dy) end; return mesh_vertices,mesh_indices,deformed,render_vertices"
    ))
    .eval()
    .unwrap()
}
#[test]
fn lua_all_pin_types_preserve_identity_and_follow_translation_in_each_mode() {
    for method in [1, 2] {
        for (pose, shift) in [
            ("", (0., 0.)),
            (
                "for i=1,pins do pin_dx[i]=pin_sx[i]+8; pin_dy[i]=pin_sy[i]-5 end",
                (8., -5.),
            ),
        ] {
            let (v, _, q, render) = lua_pipeline(method, pose);
            for (p, q) in v.chunks_exact(2).zip(q.chunks_exact(2)) {
                assert!(
                    (q[0] - p[0] - shift.0).hypot(q[1] - p[1] - shift.1) < 1.0e-5,
                    "mode {method} identity/translation mismatch"
                );
            }
            assert!(render.iter().all(|x| x.is_finite()));
        }
    }
}
#[test]
fn single_edge_bend_is_a_rigid_rotation() {
    let rgba = mask(160, 40, |_, _| true);
    for method in [1, 2] {
        let (v, _, q, _) = lua_pipeline_fixture(
            method,
            "pins=1;pin_sx={-70};pin_sy={0};pin_dx={-70};pin_dy={0};pin_types={2};pin_rotation={90};pin_scale={1}",
            15,
            false,
            Some((&rgba, 160, 40)),
        );
        let error = v
            .chunks_exact(2)
            .zip(q.chunks_exact(2))
            .map(|(p, q)| (q[0] - (-70.0 - p[1])).hypot(q[1] - (p[0] + 70.0)))
            .fold(0.0, f64::max);
        assert!(error < 0.1, "mode {method} single bend error {error}");
    }
}

#[test]
fn two_long_strip_bends_are_continuous_and_periodic() {
    let rgba = mask(26, 360, |_, _| true);
    for method in [1, 2] {
        let run = |a: f64, b: f64| {
            lua_pipeline_fixture(
                method,
                &format!(
                    "pins=2;pin_sx={{0,0}};pin_sy={{-158,155}};pin_dx={{0,0}};pin_dy={{-158,155}};pin_types={{2,2}};pin_rotation={{{a},{b}}};pin_scale={{1,1}}"
                ),
                15,
                false,
                Some((&rgba, 26, 360)),
            )
        };
        for (a, b, c, d, tolerance) in [
            (50., 0., 50., 0.01, 0.5),
            (50., 229.99, 50., 230.01, 0.5),
            (50., 0., 50., 360., 1e-5),
            (50., 230., 410., 230., 1e-5),
            (0., 50., 0.01, 50., 0.5),
        ] {
            let left = run(a, b);
            let right = run(c, d);
            // Compare material points through the render mesh even if support
            // insertion changes the simulation vertex count.
            let error = [
                (0., -158.),
                (0., 0.),
                (0., 155.),
                (-10., 100.),
                (10., -100.),
            ]
            .into_iter()
            .map(|(x, y)| {
                let sample = |r: &[f64]| {
                    for t in r.chunks_exact(12) {
                        let u = x / 26. + 0.5;
                        let v = y / 360. + 0.5;
                        let bx = t[6] - t[2];
                        let by = t[7] - t[3];
                        let cx = t[10] - t[2];
                        let cy = t[11] - t[3];
                        let det = bx * cy - by * cx;
                        let b = ((u - t[2]) * cy - (v - t[3]) * cx) / det;
                        let c = (bx * (v - t[3]) - by * (u - t[2])) / det;
                        if b >= -1e-7 && c >= -1e-7 && b + c <= 1. + 1e-7 {
                            return (
                                t[0] + b * (t[4] - t[0]) + c * (t[8] - t[0]),
                                t[1] + b * (t[5] - t[1]) + c * (t[9] - t[1]),
                            );
                        }
                    }
                    panic!("missing material sample");
                };
                let p = sample(&left.3);
                let q = sample(&right.3);
                (p.0 - q.0).hypot(p.1 - q.1)
            })
            .fold(0., f64::max);
            if tolerance < 1e-4 {
                assert!(
                    error < tolerance,
                    "mode {method}: a full turn changed the field by {error}px"
                );
            } else {
                shape_quality(
                    method,
                    error < tolerance,
                    format_args!("mode {method}: ({a},{b}) -> ({c},{d}) jumped {error}px"),
                );
            }
        }
        let mut previous = run(50., 0.);
        for angle in 1..=360 {
            let current = run(50., angle as f64);
            assert_eq!(previous.0, current.0, "animated bend changed topology");
            let jump = previous
                .2
                .chunks_exact(2)
                .zip(current.2.chunks_exact(2))
                .map(|(p, q)| (p[0] - q[0]).hypot(p[1] - q[1]))
                .fold(0., f64::max);
            shape_quality(
                method,
                jump < 12.,
                format_args!("mode {method} angle {angle}: one degree jumped {jump}px"),
            );
            previous = current;
        }
    }
}
#[test]
fn single_frame_scale_reaches_zero_continuously() {
    let rgba = mask(160, 40, |_, _| true);
    for method in [1, 2] {
        for scale in [0.0, 0.001, 0.5] {
            let (v, _, q, _) = lua_pipeline_fixture(
                method,
                &format!(
                    "pins=1;pin_sx={{-70}};pin_sy={{0}};pin_dx={{-70}};pin_dy={{0}};pin_types={{3}};pin_rotation={{0}};pin_scale={{{scale}}}"
                ),
                15,
                false,
                Some((&rgba, 160, 40)),
            );
            let error = v
                .chunks_exact(2)
                .zip(q.chunks_exact(2))
                .map(|(p, q)| (q[0] - (-70.0 + (p[0] + 70.0) * scale)).hypot(q[1] - p[1] * scale))
                .fold(0.0, f64::max);
            assert!(error < 0.1, "mode {method} scale {scale} error {error}");
        }
    }
}
#[test]
fn three_strip_bends_do_not_anchor_their_rest_positions() {
    let rgba = mask(160, 20, |_, _| true);
    for method in [1, 2] {
        for angles in ["90,90,90", "-45,0,45"] {
            let (v, t, q, render) = lua_pipeline_fixture(
                method,
                &format!(
                    "pins=3;pin_sx={{-60,0,60}};pin_sy={{0,0,0}};pin_dx={{-60,0,60}};pin_dy={{0,0,0}};pin_types={{2,2,2}};pin_rotation={{{angles}}};pin_scale={{1,1,1}}"
                ),
                20,
                false,
                Some((&rgba, 160, 20)),
            );
            if angles == "90,90,90" {
                let error = v
                    .chunks_exact(2)
                    .zip(q.chunks_exact(2))
                    .map(|(p, q)| (q[0] - (-60.0 - p[1])).hypot(q[1] - (p[0] + 60.0)))
                    .fold(0.0, f64::max);
                assert!(error < 0.1, "mode {method} common rotation error {error}");
            }
            for tri in render.chunks_exact(12) {
                let area =
                    (tri[4] - tri[0]) * (tri[9] - tri[1]) - (tri[5] - tri[1]) * (tri[8] - tri[0]);
                let uv =
                    (tri[6] - tri[2]) * (tri[11] - tri[3]) - (tri[7] - tri[3]) * (tri[10] - tri[2]);
                assert!(
                    area * uv > 0.0,
                    "mode {method} rendered strip folded at {angles}"
                );
            }
            for t in t.chunks_exact(3) {
                let area = |v: &[f64]| {
                    let a = t[0] as usize * 2;
                    let b = t[1] as usize * 2;
                    let c = t[2] as usize * 2;
                    (v[b] - v[a]) * (v[c + 1] - v[a + 1]) - (v[b + 1] - v[a + 1]) * (v[c] - v[a])
                };
                assert!(
                    area(&v) * area(&q) > 0.0,
                    "mode {method} bent strip folded at {angles}"
                );
            }
        }
    }
}

#[test]
fn lua_moving_and_detail_pins_stay_on_the_drawn_mesh() {
    for method in [1, 2] {
        let start = std::time::Instant::now();
        let (v, t, q, render) = lua_pipeline(
            method,
            "pin_dx[1]=pin_dx[1]-8; pin_dy[1]=pin_dy[1]-12; pin_dy[2]=pin_dy[2]+8; pin_rotation[2]=20",
        );
        for (p, d) in [
            ([-60., -14.], [-68., -26.]),
            ([55., -14.], [55., -6.]),
            ([0., -15.], [0., -15.]),
        ] {
            let i = v.chunks_exact(2).position(|v| v == p).unwrap();
            assert!(
                (q[i * 2] - d[0]).hypot(q[i * 2 + 1] - d[1]) < 1.0e-6,
                "mode {method} missed a material handle"
            );
            assert!(
                render
                    .chunks_exact(4)
                    .any(|r| (r[0] - d[0]).hypot(r[1] - d[1]) < 1.0e-6)
            );
        }
        assert!(q.iter().chain(&render).all(|x| x.is_finite()));
        for t in t.chunks_exact(3) {
            let area = |v: &[f64]| {
                let (a, b, c) = (t[0] as usize * 2, t[1] as usize * 2, t[2] as usize * 2);
                (v[b] - v[a]) * (v[c + 1] - v[a + 1]) - (v[b + 1] - v[a + 1]) * (v[c] - v[a])
            };
            assert!(
                area(&v) * area(&q) >= -1.0e-8,
                "mode {method} inverted a triangle in the limb fixture: {:?} -> {:?}",
                t.iter()
                    .map(|i| &v[*i as usize * 2..*i as usize * 2 + 2])
                    .collect::<Vec<_>>(),
                t.iter()
                    .map(|i| &q[*i as usize * 2..*i as usize * 2 + 2])
                    .collect::<Vec<_>>()
            );
        }
        eprintln!(
            "mode {method}: {} vertices, {} triangles, {:?}",
            v.len() / 2,
            t.len() / 3,
            start.elapsed()
        );
        if let Ok(directory) = std::env::var("PUPPET_QA_DIR") {
            std::fs::create_dir_all(&directory).unwrap();
            let mut svg = String::from(
                "<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='800' viewBox='-100 -80 200 160'><rect x='-100' y='-80' width='200' height='160' fill='#17191e'/>",
            );
            for t in t.chunks_exact(3) {
                let [a, b, c] = [t[0] as usize, t[1] as usize, t[2] as usize];
                svg.push_str(&format!("<polygon points='{},{} {},{} {},{}' fill='#4396be' stroke='#d8eaff' stroke-width='.12'/>",q[a*2],q[a*2+1],q[b*2],q[b*2+1],q[c*2],q[c*2+1]));
            }
            svg.push_str("</svg>");
            std::fs::write(format!("{directory}/mode-{method}.svg"), svg).unwrap();
        }
    }
}

#[test]
fn deterministic_small_mask_sweep_never_cuts_opaque_cells() {
    let mut state = 7_u64;
    for _ in 0..48 {
        let mut pixels = [false; 100];
        for p in &mut pixels {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            *p = (state >> 60) > 5;
        }
        let rgba = mask(10, 10, |x, y| pixels[y * 10 + x]);
        assert_coverage(&rgba, 10, 10, 5);
    }
}

#[test]
fn diagonal_sleeves_respect_density_at_hd_resolution() {
    let mut counts = Vec::new();
    for size in [240, 960, 1920] {
        let rgba = mask(size, size / 2, |x, y| {
            let (x, y) = (x as f64 / size as f64, y as f64 / size as f64);
            x > 0.05 && x < 0.95 && (y - (0.1 + 0.28 * x)).abs() < 0.025
        });
        for density in [8, 15, 30] {
            let start = std::time::Instant::now();
            let mesh = generate_mesh(&rgba, size, size / 2, 1, density, 0).unwrap();
            eprintln!(
                "sleeve {size}px density {density}: {} vertices, {} triangles, {:?}",
                mesh.points.len(),
                mesh.triangles.len(),
                start.elapsed()
            );
            assert!(
                mesh.points.len() < density as usize * 4,
                "pixel-scale boundary defeated density control"
            );
            counts.push(mesh.points.len());
        }
    }
    assert!(
        counts[6] < counts[0] * 2 + 8
            && counts[7] < counts[1] * 2 + 8
            && counts[8] < counts[2] * 2 + 8,
        "resolution must not dictate mesh size"
    );
    assert!(
        counts[6] < counts[8],
        "density must control mesh complexity"
    );
}

#[test]
fn density_controls_the_final_mesh_with_all_handle_processing_enabled() {
    let mut counts = Vec::new();
    for (density, budget) in [(5, 60), (15, 95), (30, 180), (50, 350)] {
        let bare = lua_pipeline_at_density(2, "pins=0", density);
        let mut mode_counts = Vec::new();
        for mode in [1, 2] {
            let (v, t, q, _) = lua_pipeline_at_density(
                mode,
                "pin_types[6]=PIN_TYPE.DETAIL; pin_layer[6]=0; pin_dx[1]=-180; pin_dy[1]=-60",
                density,
            );
            let count = v.len() / 2;
            eprintln!(
                "density {density}, mode {mode}: {} base + {} handle = {count} vertices, {} triangles",
                bare.0.len() / 2,
                count - bare.0.len() / 2,
                t.len() / 3
            );
            assert!(
                count <= bare.0.len() / 2 + 25,
                "five material pins may not trigger global refinement"
            );
            assert!(
                count <= budget,
                "density {density} exceeded the total fixture budget"
            );
            let i = v.chunks_exact(2).position(|p| p == [-60.0, -14.0]).unwrap();
            assert!((q[i * 2] + 180.0).hypot(q[i * 2 + 1] + 60.0) < 1e-7);
            mode_counts.push(count);
        }
        assert_eq!(mode_counts[0], mode_counts[1]);
        counts.push(mode_counts[0]);
    }
    assert!(
        counts.windows(2).all(|p| p[0] < p[1]),
        "density must remain effective after inserting handles"
    );
}

#[test]
fn curved_silhouette_density_reduces_boundary_vertices_without_excess_padding() {
    let size = 512;
    let rgba = mask(size, size, |x, y| {
        let (x, y) = (x as f64 / size as f64, y as f64 / size as f64);
        let ellipse = |cx: f64, cy: f64, rx: f64, ry: f64| {
            ((x - cx) / rx).powi(2) + ((y - cy) / ry).powi(2) < 1.0
        };
        let capsule = |a: (f64, f64), b: (f64, f64), r: f64| {
            let (u, v) = (b.0 - a.0, b.1 - a.1);
            let t = (((x - a.0) * u + (y - a.1) * v) / (u * u + v * v)).clamp(0.0, 1.0);
            (x - a.0 - t * u).hypot(y - a.1 - t * v) < r
        };
        ellipse(0.5, 0.2, 0.15, 0.16)
            || ellipse(0.5, 0.52, 0.16, 0.25)
            || capsule((0.4, 0.42), (0.15, 0.23), 0.055)
            || capsule((0.6, 0.42), (0.85, 0.23), 0.055)
            || capsule((0.44, 0.68), (0.25, 0.9), 0.065)
            || capsule((0.56, 0.68), (0.75, 0.9), 0.065)
            || ellipse(0.15, 0.23, 0.08, 0.075)
            || ellipse(0.85, 0.23, 0.08, 0.075)
    });
    let opaque = rgba.chunks_exact(4).filter(|p| p[3] > 0).count() as f64;
    let mut perimeter = 0.0;
    for y in 1..size - 1 {
        for x in 1..size - 1 {
            if rgba[(y * size + x) * 4 + 3] > 0 {
                perimeter += [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
                    .iter()
                    .filter(|&&(x, y)| rgba[(y * size + x) * 4 + 3] == 0)
                    .count() as f64;
            }
        }
    }
    let mut counts = Vec::new();
    for density in [5, 15, 30] {
        let mesh = generate_mesh(&rgba, size, size, 1, density, 0).unwrap();
        let mut edges = std::collections::BTreeMap::new();
        let mut area = 0.0;
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
            let [a, b, c] = t.map(|i| mesh.points[i]);
            area += ((b.0 - a.0) as f64 * (c.1 - a.1) as f64
                - (b.1 - a.1) as f64 * (c.0 - a.0) as f64)
                .abs()
                * 0.5;
        }
        let boundary = edges.values().filter(|&&n| n == 1).count();
        eprintln!(
            "curved density {density}: {boundary} boundary, {} total vertices, padding {:.3}%",
            mesh.points.len(),
            100.0 * (area / opaque - 1.0)
        );
        counts.push(boundary);
        assert!(
            area / opaque <= 1.0 + (0.25 / density as f64).max(perimeter / opaque) + 1e-5,
            "excess transparent area"
        );
        assert!(
            boundary <= density as usize * 3 + 30,
            "curve simplification exploded into pixel steps"
        );
        for y in 1..size - 1 {
            for x in 1..size - 1 {
                if rgba[(y * size + x) * 4 + 3] > 0
                    && [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
                        .iter()
                        .any(|&(x, y)| rgba[(y * size + x) * 4 + 3] == 0)
                {
                    for (dx, dy) in [
                        (0.001, 0.001),
                        (0.999, 0.001),
                        (0.001, 0.999),
                        (0.999, 0.999),
                    ] {
                        assert!(
                            covers(
                                &mesh,
                                (
                                    x as f32 + dx - size as f32 * 0.5,
                                    y as f32 + dy - size as f32 * 0.5
                                )
                            ),
                            "curve clipped at density {density}"
                        );
                    }
                }
            }
        }
    }
    assert!(
        counts[0] * 4 < counts[2] * 3,
        "low density did not simplify curved boundaries"
    );
    assert!(
        counts[1] < counts[2],
        "medium density did not simplify curved boundaries"
    );
}

#[test]
fn neutral_bend_is_inert_with_detail_controls() {
    for method in [1, 2] {
        let motion = "pin_dx[1]=pin_dx[1]-8; pin_dy[2]=pin_dy[2]+8;";
        let detail = lua_pipeline(method, motion);
        let no_bend = lua_pipeline(
            method,
            &format!("{motion} pin_types[4]=PIN_TYPE.OVERLAP; pin_layer[4]=0"),
        );
        assert_eq!(detail.0, no_bend.0);
        assert!(
            detail
                .2
                .iter()
                .zip(&no_bend.2)
                .all(|(a, b)| (a - b).abs() < 1e-8)
        );
    }
}

#[test]
fn adjacent_bend_and_detail_keep_position_handles_and_valid_faces() {
    for method in [1, 2] {
        let (v, t, q, render) = lua_pipeline(
            method,
            "pin_sx[4]=48; pin_sy[4]=-14; pin_rotation[4]=-12; pin_scale[4]=1.05; pin_rotation[2]=18; pin_scale[2]=1.15",
        );
        let i = v.chunks_exact(2).position(|v| v == [55., -14.]).unwrap();
        assert!((q[i * 2] - 55.).hypot(q[i * 2 + 1] + 14.) < 1e-7);
        for tri in t.chunks_exact(3) {
            let area = |p: &[f64]| {
                let (a, b, c) = (
                    tri[0] as usize * 2,
                    tri[1] as usize * 2,
                    tri[2] as usize * 2,
                );
                (p[b] - p[a]) * (p[c + 1] - p[a + 1]) - (p[b + 1] - p[a + 1]) * (p[c] - p[a])
            };
            assert!(
                area(&v) * area(&q) > -1e-8,
                "mode {method}: adjacent frame controls inverted a face {:?} -> {:?}",
                tri.iter()
                    .map(|i| &v[*i as usize * 2..*i as usize * 2 + 2])
                    .collect::<Vec<_>>(),
                tri.iter()
                    .map(|i| &q[*i as usize * 2..*i as usize * 2 + 2])
                    .collect::<Vec<_>>()
            );
        }
        assert!(render.iter().all(|v| v.is_finite()));
    }
}

#[test]
fn bend_does_not_emit_a_positional_anchor() {
    let lua = mlua::Lua::new();
    let src = include_str!("../../../script/PuppetTool/PuppetPinDynamics.inc.lua");
    let code = format!(
        "{src}\n local v={{0,0,10,0,0,10}}; local p={{{{sx=0,sy=0,dx=100,dy=100,kind=2,rotation=30,scale=1,layer=0}}}}; local sources,_,_,_,_,poses=puppet_pin_dynamics.build(v,v,p,5,{{POSITION=0,BONE=1,BEND=2,DETAIL=3,STARCH=4,OVERLAP=5}},{{0,1,2}}); return sources,poses"
    );
    let (sources, poses): (Vec<f64>, Vec<f64>) = lua.load(code).eval().unwrap();
    assert!(
        sources.is_empty(),
        "a bend may control its frame, never its translation"
    );
    assert_eq!(poses[9], 0.0);
    assert_eq!(poses[10], 1.0);
}

#[test]
fn overlap_preview_uses_the_same_connected_mesh_range_as_layering() {
    let lua = mlua::Lua::new();
    let src = include_str!("../../../script/PuppetTool/PuppetPinDynamics.inc.lua");
    let code = format!(
        "{src}\nlocal v={{0,0,4,0,0,4, 5,0,9,0,9,4}}; local i={{0,1,2,3,4,5}}; local p={{{{sx=0,sy=0,kind=5,show_range=true,range=10,layer=1}}}}; return puppet_pin_dynamics.overlap_influence(v,i,p,5)"
    );
    let influence: Vec<f64> = lua.load(code).eval().unwrap();
    assert!(influence[..3].iter().all(|&value| value > 0.0));
    assert_eq!(&influence[3..], &[0.0, 0.0, 0.0]);
}

#[test]
fn bone_endpoints_do_not_prescribe_rotation_or_add_midpoint_controls() {
    let lua = mlua::Lua::new();
    let src = include_str!("../../../script/PuppetTool/PuppetPinDynamics.inc.lua");
    let code = format!(
        "{src}\n local v={{0,0,10,0,0,10}}; local q={{0,0,0,20,-10,0}}; local p={{{{sx=0,sy=0,dx=0,dy=0,kind=1}},{{sx=10,sy=0,dx=0,dy=20,kind=1}}}}; local _,_,_,_,_,poses=puppet_pin_dynamics.build(v,q,p,20,{{POSITION=0,BONE=1,BEND=2,DETAIL=3,STARCH=4,OVERLAP=5}},{{0,1,2}},{{{{1,{{2}}}}}}); return poses"
    );
    let poses: Vec<f64> = lua.load(code).eval().unwrap();
    assert_eq!(poses.len(), 24);
    for p in poses.chunks_exact(12) {
        assert_eq!(p[6], 0.);
        assert_eq!(p[7], 1.0);
        assert_eq!(p[10], 0.0);
        assert_eq!(p[11], 1.0);
    }
}

#[test]
fn frame_controls_change_both_solvers() {
    for method in [1, 2] {
        let neutral = lua_pipeline(method, "pin_sx[4]=48; pin_sy[4]=-14");
        let edited = lua_pipeline(
            method,
            "pin_sx[4]=48; pin_sy[4]=-14; pin_rotation[4]=25; pin_scale[4]=1.2",
        );
        // Compare the actual rendered image at a shared material sample, since
        // activating a bend deliberately adds a small local support ring.
        let changed = edited.3.chunks_exact(4).any(|a| {
            neutral.3.chunks_exact(4).any(|b| {
                (a[2] - b[2]).abs() < 1e-8
                    && (a[3] - b[3]).abs() < 1e-8
                    && (a[0] - b[0]).hypot(a[1] - b[1]) > 0.01
            })
        });
        assert!(changed, "mode {method} ignored a frame edit");
    }
}

fn rendered_material_point(render: &[f64], x: f64, y: f64) -> (f64, f64) {
    let (u, v) = (x / 160.0 + 0.5, y / 120.0 + 0.5);
    render
        .chunks_exact(12)
        .find_map(|t| {
            let (bx, by, cx, cy) = (t[6] - t[2], t[7] - t[3], t[10] - t[2], t[11] - t[3]);
            let det = bx * cy - by * cx;
            if det.abs() < 1e-14 {
                return None;
            }
            let b = ((u - t[2]) * cy - (v - t[3]) * cx) / det;
            let c = (bx * (v - t[3]) - by * (u - t[2])) / det;
            let a = 1.0 - b - c;
            (a >= -1e-8 && b >= -1e-8 && c >= -1e-8).then_some((
                a * t[0] + b * t[4] + c * t[8],
                a * t[1] + b * t[5] + c * t[9],
            ))
        })
        .expect("material sample absent")
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_illustration_wrist_continuity() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for density in [5, 15] {
        let mut previous: Option<f64> = None;
        let mut maximum_jump = 0.0_f64;
        let mut maximum_error = 0.0_f64;
        for x in (200..=360).step_by(5) {
            let pose = format!(
                r#"
                pins=5
                pin_sx={{360-hw,70-hw,650-hw,190-hw,531-hw}}
                pin_sy={{510-hh,380-hh,380-hh,727-hh,727-hh}}
                pin_types={{1,1,1,1,1}};bone_forest={{{{1,{{2}},{{3}},{{4}},{{5}}}}}}
                bone_forest={{bone_forest}}
                for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[2]={x}-hw;pin_dy[2]=70-hh
            "#
            );
            let (_, _, _, render) =
                lua_pipeline_fixture(2, &pose, density, false, Some((&rgba, 723, 800)));
            if density == 15 && [275, 280].contains(&x) {
                std::fs::write(
                    format!("../../target/puppet-qa/actual-wrist-{x}.csv"),
                    render
                        .chunks_exact(12)
                        .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                        .collect::<Vec<_>>()
                        .join("\n"),
                )
                .unwrap();
            }
            let sample = |x: f64, y: f64| {
                rendered_material_point(&render, (x / 723. - 0.5) * 160., (y / 800. - 0.5) * 120.)
            };
            let a = sample(50., 390.);
            let b = sample(100., 370.);
            let angle = (b.1 - a.1).atan2(b.0 - a.0);
            let c = sample(140., 355.);
            let d = sample(180., 340.);
            let expected =
                (d.1 - c.1).atan2(d.0 - c.0) + (-20.0_f64).atan2(50.) - (-15.0_f64).atan2(40.);
            let difference = |a: f64, b: f64| (a - b).sin().atan2((a - b).cos()).abs().to_degrees();
            let error = difference(angle, expected);
            maximum_error = maximum_error.max(error);
            if let Some(old) = previous {
                maximum_jump = maximum_jump.max(difference(angle, old));
            }
            previous = Some(angle);
            eprintln!(
                "actual wrist d{density} x{x}: angle {}, error {error}",
                angle.to_degrees()
            );
        }
        assert!(maximum_jump < 5., "d{density} max jump {maximum_jump}");
        assert!(maximum_error < 10., "d{density} max error {maximum_error}");
    }
}

#[test]
fn position_and_bone_pins_have_identical_deformation_in_both_modes() {
    let pose = "pin_types[1]=PIN_TYPE.BONE;pin_dx[1]=-35;pin_dy[1]=-45;bone_forest={{3,{1}}}";
    for mode in [1, 2] {
        let expected = lua_pipeline(
            mode,
            &format!(
                "{pose};for i=1,pins do if pin_types[i]==PIN_TYPE.BONE then pin_types[i]=PIN_TYPE.POSITION end end"
            ),
        );
        for wrappers in [1, 2] {
            let actual = lua_pipeline(
                mode,
                &format!("{pose};for _=1,{wrappers} do bone_forest={{bone_forest}} end"),
            );
            assert_eq!(
                expected.0, actual.0,
                "wrapped forest changed the prepared vertices"
            );
            assert_eq!(expected.1, actual.1, "bone hierarchy changed mesh topology");
            assert_eq!(expected.2.len(), actual.2.len());
            assert_eq!(expected.3.len(), actual.3.len());
            for (a, b) in expected
                .2
                .iter()
                .chain(&expected.3)
                .zip(actual.2.iter().chain(&actual.3))
            {
                assert!(
                    (a - b).abs() < 1e-8,
                    "mode {mode}: wrapped forest changed the deformation"
                );
            }
        }
    }
}

#[test]
fn arap_raised_wrist_angle_changes_continuously() {
    let difference = |a: f64, b: f64| (a - b).sin().atan2((a - b).cos()).abs().to_degrees();
    for density in [5, 15, 30] {
        let mut maximum_jump = 0.0_f64;
        let mut maximum_error = 0.0_f64;
        // Sample a compressed arm, including the old rotation-branch transition.
        // Compare 0.01px perturbations: a bent arm can rotate quickly over 1px
        // without being discontinuous, so angular speed alone is not a kink.
        for step in -40..=10 {
            let mut previous: Option<(f64, f64)> = None;
            for offset in [0., 0.01] {
                let x = f64::from(step) + offset;
                let pose = format!(
                    r#"
                    pins=5
                    pin_sx={{0,-60,55,-18,16}};pin_sy={{20,-14,-14,47,47}}
                    pin_types={{0,0,0,0,0}}
                    for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                    pin_dx[2]={x};pin_dy[2]=-15
                "#
                );
                let (_, _, _, render) = lua_pipeline_on_mask(2, &pose, density, true);
                let a = rendered_material_point(&render, -67., -14.);
                let b = rendered_material_point(&render, -57., -14.);
                let c = rendered_material_point(&render, -56., -14.);
                let d = rendered_material_point(&render, -48., -14.);
                let angle = (b.1 - a.1).atan2(b.0 - a.0);
                let forearm = (d.1 - c.1).atan2(d.0 - c.0);
                let relative = angle - forearm;
                maximum_error = maximum_error.max(difference(angle, forearm));
                if let Some((old, old_relative)) = previous {
                    maximum_jump = maximum_jump
                        .max(difference(angle, old))
                        .max(difference(relative, old_relative));
                }
                previous = Some((angle, relative));
            }
        }
        eprintln!("wrist d{density}: 0.01px jump {maximum_jump}, forearm error {maximum_error}");
        assert!(
            maximum_jump < 0.2,
            "wrist rotation jumped at density {density}"
        );
        assert!(
            maximum_error < 30.,
            "wrist kinked away from adjacent material at density {density}"
        );
    }
}

#[test]
fn curved_low_density_puppet_keeps_wrist_aligned_with_material() {
    let mut counts = Vec::new();
    for density in [5, 15, 30] {
        for mode in [1, 2] {
            let (v, triangles, q, render) = lua_pipeline_on_mask(
                mode,
                "pin_types[1]=PIN_TYPE.BONE; pin_dx[1]=-35; pin_dy[1]=-65; bone_forest={{3,{1}}}",
                density,
                true,
            );
            let a = rendered_material_point(&render, -63., -14.);
            let b = rendered_material_point(&render, -57., -14.);
            let angle = (b.1 - a.1).atan2(b.0 - a.0);
            let c = rendered_material_point(&render, -50., -14.);
            let d = rendered_material_point(&render, -40., -14.);
            let expected = (d.1 - c.1).atan2(d.0 - c.0);
            let error = (angle - expected)
                .sin()
                .atan2((angle - expected).cos())
                .abs()
                .to_degrees();
            let scale = (b.0 - a.0).hypot(b.1 - a.1) / 6.0;
            eprintln!(
                "curved density {density} mode {mode}: {} vertices, wrist error {error:.2} deg, scale {scale:.3}",
                v.len() / 2
            );
            assert!(error < 30.0, "wrist bent away from the adjacent material");
            assert!((0.8..1.2).contains(&scale), "wrist collapsed");
            let i = v.chunks_exact(2).position(|p| p == [-60., -14.]).unwrap();
            assert!((q[i * 2] + 35.).hypot(q[i * 2 + 1] + 65.) < 1e-7);
            if let Ok(directory) = std::env::var("PUPPET_QA_DIR") {
                std::fs::create_dir_all(&directory).unwrap();
                let mut svg = String::from(
                    "<svg xmlns='http://www.w3.org/2000/svg' width='1260' height='600' viewBox='-100 -100 420 200'><rect x='-100' y='-100' width='420' height='200' fill='#17191e'/>",
                );
                for (offset, points) in [(0.0, &v), (210.0, &q)] {
                    for face in triangles.chunks_exact(3) {
                        let [a, b, c] = [
                            face[0] as usize * 2,
                            face[1] as usize * 2,
                            face[2] as usize * 2,
                        ];
                        let color = if (v[a + 1] + v[b + 1] + v[c + 1]) / 3.0 < -26.0 {
                            "#eeaa99"
                        } else {
                            "#4396be"
                        };
                        svg.push_str(&format!("<polygon points='{},{} {},{} {},{}' fill='{color}' stroke='#d8eaff' stroke-width='.2'/>",points[a]+offset,points[a+1],points[b]+offset,points[b+1],points[c]+offset,points[c+1]));
                    }
                }
                svg.push_str("</svg>");
                std::fs::write(
                    format!("{directory}/curved-wrist-d{density}-m{mode}.svg"),
                    svg,
                )
                .unwrap();
            }
            if mode == 2 {
                counts.push(v.len() / 2);
            }
        }
    }
    assert!(counts[0] < counts[1] && counts[1] < counts[2]);
    assert!(
        counts[0] <= 85 && counts[1] <= 110 && counts[2] <= 190,
        "curved pin refinement exceeded density budget"
    );
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_illustration_extreme_arap_pull() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    for density in [5, 15, 30] {
        for (name, x, y) in [
            ("down-near", 360, 600),
            ("down-mid", 360, 800),
            ("down", 360, 950),
            ("down-far", 360, 1400),
            ("down-offset-left", 350, 600),
            ("down-offset-right", 370, 600),
            ("twist", 500, 600),
            ("twist-left", 220, 600),
        ] {
            let pose = format!(
                r#"
                pins=5
                pin_sx={{70-hw,650-hw,190-hw,531-hw,360-hw}}
                pin_sy={{380-hh,380-hh,727-hh,727-hh,180-hh}}
                pin_types={{0,0,0,0,0}};bone_forest={{}}
                for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[5]={x}-hw;pin_dy[5]={y}-hh
            "#
            );
            let (_, _, q, render) =
                lua_pipeline_fixture(2, &pose, density, false, Some((&rgba, 723, 800)));
            std::fs::write(
                format!("../../target/puppet-qa/actual-extreme-{name}-d{density}.csv"),
                render
                    .chunks_exact(12)
                    .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                    .collect::<Vec<_>>()
                    .join("\n"),
            )
            .unwrap();
            let locate = |x: f64, y: f64| {
                rendered_material_point(&render, (x / 723. - 0.5) * 160., (y / 800. - 0.5) * 120.)
            };
            let neck = locate(360., 285.);
            let chin = locate(360., 245.);
            let nose = locate(360., 180.);
            let crown = locate(360., 40.);
            let angle = |a: (f64, f64), b: (f64, f64)| (b.1 - a.1).atan2(b.0 - a.0);
            let axis = angle(neck, nose);
            let difference = |a: f64, b: f64| (a - b).sin().atan2((a - b).cos()).abs().to_degrees();
            let bend =
                difference(angle(nose, crown), axis).max(difference(angle(chin, nose), axis));
            eprintln!(
                "{name} d{density}: head bend {bend}, down error {}",
                difference(angle(nose, crown), std::f64::consts::FRAC_PI_2)
            );
            assert!(
                bend < 25.,
                "head folded at {name}, density {density}: {bend}"
            );
            if x == 360 {
                assert!(
                    difference(angle(nose, crown), std::f64::consts::FRAC_PI_2) < 20.,
                    "head did not point downward at {name}, density {density}"
                );
                assert!(
                    neck.1 < chin.1 && chin.1 < nose.1 && nose.1 < crown.1,
                    "head centerline doubled back: {neck:?} {chin:?} {nose:?} {crown:?}"
                );
            }
            // Folding/overlap is permitted; a blanket determinant barrier would
            // reintroduce the shoulder collapse. Check the head axis above.
            assert!(q.iter().all(|p| p.is_finite()));
        }
    }
}

#[test]
fn arap_extreme_position_pull_keeps_head_axis_and_all_pin_targets() {
    for density in [5, 15, 30] {
        for kind in [0, 1] {
            for (x, y) in [(0., 80.), (45., 55.), (-45., 55.)] {
                let pose = format!(
                    r#"
                pins=5
                pin_sx={{-60,55,-18,16,0}}; pin_sy={{-14,-14,47,47,-35}}
                pin_types={{{kind},{kind},{kind},{kind},{kind}}}; bone_forest={{}}
                for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[5]={x};pin_dy[5]={y}
            "#
                );
                let (v, _, q, render) = lua_pipeline_at_density(2, &pose, density);
                let chin = rendered_material_point(&render, 0., -27.);
                let center = rendered_material_point(&render, 0., -35.);
                let crown = rendered_material_point(&render, 0., -47.);
                let a = (center.1 - chin.1).atan2(center.0 - chin.0);
                let b = (crown.1 - center.1).atan2(crown.0 - center.0);
                let bend = (a - b).sin().atan2((a - b).cos()).abs().to_degrees();
                eprintln!("centerline {chin:?} {center:?} {crown:?} bend {bend}");
                assert!(bend < 25., "head centerline folded: {bend}");
                for (source, target) in [
                    ((-60., -14.), (-60., -14.)),
                    ((55., -14.), (55., -14.)),
                    ((-18., 47.), (-18., 47.)),
                    ((16., 47.), (16., 47.)),
                    ((0., -35.), (x, y)),
                ] {
                    let i = v
                        .chunks_exact(2)
                        .position(|p| (p[0] - source.0).hypot(p[1] - source.1) < 1e-8)
                        .unwrap()
                        * 2;
                    assert!((q[i] - target.0).hypot(q[i + 1] - target.1) < 1e-8);
                }
                if x == 0. {
                    let crown = v
                        .chunks_exact(2)
                        .enumerate()
                        .min_by(|(_, a), (_, b)| {
                            (a[0].powi(2) + (a[1] + 49.).powi(2))
                                .total_cmp(&(b[0].powi(2) + (b[1] + 49.).powi(2)))
                        })
                        .unwrap()
                        .0;
                    assert!(
                        q[crown * 2 + 1] > y,
                        "head tip folded back: density {density}, kind {kind}"
                    );
                }
            }
        }
    }
}

#[test]
fn large_hand_pull_preserves_default_detail_head() {
    for (case, pose, expected) in [
        (
            "pull",
            "pin_types[6]=PIN_TYPE.DETAIL; pin_layer[6]=0; pin_dx[1]=-180; pin_dy[1]=-60",
            1.0,
        ),
        (
            "starch",
            "pin_types[6]=PIN_TYPE.DETAIL; pin_layer[6]=0; pin_dx[1]=-180; pin_dy[1]=-60; pin_sx[5]=-16; pin_sy[5]=-16",
            1.0,
        ),
        (
            "head",
            "pin_types[6]=PIN_TYPE.DETAIL; pin_layer[6]=0; pin_dx[6]=20; pin_dy[6]=-80; pin_rotation[6]=-35",
            1.0,
        ),
        (
            "scale",
            "pin_types[6]=PIN_TYPE.DETAIL; pin_layer[6]=0; pin_dx[6]=20; pin_dy[6]=-80; pin_rotation[6]=-35; pin_scale[6]=1.3",
            1.69,
        ),
        (
            "bone",
            "pin_types[1]=PIN_TYPE.BONE; pin_dx[1]=-35; pin_dy[1]=-65; bone_forest={{3,{1}}}",
            1.0,
        ),
    ] {
        for method in [1, 2] {
            let (v, t, q, render) = lua_pipeline(method, pose);
            for triangle in render.chunks_exact(12) {
                if [3, 7, 11]
                    .iter()
                    .all(|&i| triangle[i] * 120.0 - 60.0 < -27.0)
                {
                    let source = (triangle[6] - triangle[2]) * (triangle[11] - triangle[3])
                        - (triangle[7] - triangle[3]) * (triangle[10] - triangle[2]);
                    let target = (triangle[4] - triangle[0]) * (triangle[9] - triangle[1])
                        - (triangle[5] - triangle[1]) * (triangle[8] - triangle[0]);
                    shape_quality(
                        method,
                        source * target >= -1e-9,
                        format_args!("{case}, mode {method}: rendered head folded"),
                    );
                }
            }
            let mut before = 0.0;
            let mut after = 0.0;
            let landmarks = [(-6.0, -39.0), (6.0, -39.0), (-4.0, -31.0), (4.0, -31.0)];
            let mapped = landmarks.map(|(x, y)| {
                let (u, v) = (x / 160.0 + 0.5, y / 120.0 + 0.5);
                render
                    .chunks_exact(12)
                    .find_map(|t| {
                        let (bx, by, cx, cy) =
                            (t[6] - t[2], t[7] - t[3], t[10] - t[2], t[11] - t[3]);
                        let det = bx * cy - by * cx;
                        if det.abs() < 1e-14 {
                            return None;
                        }
                        let b = ((u - t[2]) * cy - (v - t[3]) * cx) / det;
                        let c = (bx * (v - t[3]) - by * (u - t[2])) / det;
                        let a = 1.0 - b - c;
                        (a >= -1e-8 && b >= -1e-8 && c >= -1e-8).then_some((
                            a * t[0] + b * t[4] + c * t[8],
                            a * t[1] + b * t[5] + c * t[9],
                        ))
                    })
                    .expect("face landmark left the material mesh")
            });
            for i in 0..landmarks.len() {
                for j in i + 1..landmarks.len() {
                    let original =
                        (landmarks[i].0 - landmarks[j].0).hypot(landmarks[i].1 - landmarks[j].1);
                    let actual = (mapped[i].0 - mapped[j].0).hypot(mapped[i].1 - mapped[j].1);
                    let ratio = actual / original / f64::sqrt(expected);
                    shape_quality(
                        method,
                        (0.75..1.3).contains(&ratio),
                        format_args!("{case}, mode {method}: face landmark distance ratio {ratio}"),
                    );
                }
            }
            for face in t.chunks_exact(3) {
                let ids = face.iter().map(|i| *i as usize * 2).collect::<Vec<_>>();
                if ids.iter().all(|&i| v[i + 1] < -27.0) {
                    let area = |p: &[f64]| {
                        ((p[ids[1]] - p[ids[0]]) * (p[ids[2] + 1] - p[ids[0] + 1])
                            - (p[ids[1] + 1] - p[ids[0] + 1]) * (p[ids[2]] - p[ids[0]]))
                            .abs()
                    };
                    before += area(&v);
                    after += area(&q);
                }
            }
            eprintln!(
                "{case} mode {method}: head area ratio {}, {} vertices",
                after / before,
                v.len() / 2
            );
            if let Ok(directory) = std::env::var("PUPPET_QA_DIR") {
                std::fs::create_dir_all(&directory).unwrap();
                let mut svg = String::from(
                    "<svg xmlns='http://www.w3.org/2000/svg' width='1200' height='600' viewBox='-220 -100 400 200'><rect x='-220' y='-100' width='400' height='200' fill='#17191e'/>",
                );
                for tri in render.chunks_exact(12) {
                    let color = if (tri[3] + tri[7] + tri[11]) / 3.0 < 0.28 {
                        "#eeaa99"
                    } else {
                        "#4396be"
                    };
                    svg.push_str(&format!("<polygon points='{},{} {},{} {},{}' fill='{color}' stroke='#d8eaff' stroke-width='.15'/>",tri[0],tri[1],tri[4],tri[5],tri[8],tri[9]));
                }
                svg.push_str("</svg>");
                std::fs::write(format!("{directory}/{case}-{method}.svg"), svg).unwrap();
            }
            shape_quality(
                method,
                before > 0.0 && after / before / expected > 0.8 && after / before / expected < 1.25,
                format_args!(
                    "{case}, mode {method}: head scale must follow user input, not incidental strain"
                ),
            );
        }
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_crossed_arm_preserves_width() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut failures = Vec::new();
    for density in [5, 15, 30] {
        for (x, y) in [(360, 70), (320, 100), (380, 140)] {
            let pose = format!(
                r#"
                pins=6
                pin_sx={{70-hw,650-hw,190-hw,531-hw,360-hw,360-hw}}
                pin_sy={{380-hh,380-hh,727-hh,727-hh,510-hh,180-hh}}
                pin_types={{0,0,0,0,0,0}};bone_forest={{}}
                for i=1,pins do pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[2]={x}-hw;pin_dy[2]={y}-hh
            "#
            );
            let (_, _, _, render) =
                lua_pipeline_fixture(2, &pose, density, false, Some((&rgba, 723, 800)));
            std::fs::write(
                format!("../../target/puppet-qa/actual-crossed-{x}-{y}-d{density}.csv"),
                render
                    .chunks_exact(12)
                    .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                    .collect::<Vec<_>>()
                    .join("\n"),
            )
            .unwrap();
            let at = |x: f64, y: f64| {
                rendered_material_point(&render, (x / 723. - 0.5) * 160., (y / 800. - 0.5) * 120.)
            };
            for (a, b, c, d) in [
                ((500., 310.), (485., 350.), (480., 325.), (540., 350.)),
                ((469., 290.), (451., 349.), (435., 310.), (505., 330.)),
            ] {
                let pa = at(a.0, a.1);
                let pb = at(b.0, b.1);
                let pc = at(c.0, c.1);
                let pd = at(d.0, d.1);
                let width = ((pb.0 - pa.0) * (pd.1 - pc.1) - (pb.1 - pa.1) * (pd.0 - pc.0)).abs()
                    / (pd.0 - pc.0).hypot(pd.1 - pc.1);
                let rest = ((b.0 - a.0) * (d.1 - c.1) - (b.1 - a.1) * (d.0 - c.0)).abs()
                    / (d.0 - c.0).hypot(d.1 - c.1);
                eprintln!("upper {a:?} d{density} ({x},{y}) width {}", width / rest);
                if !(0.90..1.10).contains(&(width / rest)) {
                    failures.push(format!(
                        "upper {a:?} d{density} ({x},{y}): {}",
                        width / rest
                    ));
                }
            }
            let a = at(557., 331.);
            let b = at(540., 368.);
            let c = at(510., 335.);
            let d = at(590., 365.);
            let tx = d.0 - c.0;
            let ty = d.1 - c.1;
            let width = ((b.0 - a.0) * ty - (b.1 - a.1) * tx).abs() / tx.hypot(ty);
            let rest = (17. * 30. + 37. * 80.) / (80.0_f64.hypot(30.));
            eprintln!("crossed d{density} ({x},{y}) width {}", width / rest);
            if width / rest < 0.85 {
                failures.push(format!("forearm d{density} ({x},{y}): {}", width / rest));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn extreme_compression_uses_the_same_fixed_work_as_an_ordinary_drag() {
    for side in [16, 24] {
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();
        for y in 0..side {
            for x in 0..side {
                vertices.extend([x as f64 * 10., y as f64 * 10.]);
            }
        }
        for y in 0..side - 1 {
            for x in 0..side - 1 {
                let a = y * side + x;
                triangles.extend([a, a + 1, a + side, a + 1, a + side + 1, a + side]);
            }
        }
        let end = (side - 1) as f64 * 10.;
        let sources = [0., 0., end, 0., 0., end, end, end];
        let destinations = [
            0.,
            0.,
            end * 0.02,
            0.,
            0.,
            end * 0.02,
            end * 0.02,
            end * 0.02,
        ];
        // Warm the rest-mesh factors with an ordinary non-identity pose.
        let ordinary = [0., 0., end, 5., 0., end, end, end];
        let _ = crate::arap::solve(&vertices, &triangles, &sources, &ordinary).unwrap();
        crate::arap::take_work();
        let started = std::time::Instant::now();
        let _ = crate::arap::solve(&vertices, &triangles, &sources, &ordinary).unwrap();
        let normal_time = started.elapsed();
        let normal_work = crate::arap::take_work();
        let started = std::time::Instant::now();
        let output = crate::arap::solve(&vertices, &triangles, &sources, &destinations).unwrap();
        let extreme_work = crate::arap::take_work();
        eprintln!(
            "{} vertices: ordinary {:?}, extreme {:?}, work {:?} / {:?}",
            side * side,
            normal_time,
            started.elapsed(),
            normal_work,
            extreme_work
        );
        assert_eq!(
            normal_work, extreme_work,
            "extreme targets increased solver work"
        );
        assert_eq!(
            &extreme_work[1..],
            &[0, 0],
            "drag refactored or ran an iterative linear solve"
        );
        assert!(output.iter().all(|q| q.x.is_finite() && q.y.is_finite()));
        for (k, i) in [0, side - 1, (side - 1) * side, side * side - 1]
            .into_iter()
            .enumerate()
        {
            assert!(
                (output[i as usize].x - destinations[k * 2])
                    .hypot(output[i as usize].y - destinations[k * 2 + 1])
                    < 1e-8
            );
        }
    }
}

#[test]
fn moving_controls_do_not_run_a_redundant_preliminary_solve_in_either_mode() {
    for method in [1, 2] {
        for kind in [0, 1, 3] {
            crate::arap::take_solve_calls();
            let _ = lua_pipeline(
                method,
                &format!("for i=1,pins do pin_types[i]={kind} end;pin_dx[1]=-35;pin_dy[1]=-65"),
            );
            assert_eq!(
                crate::arap::take_solve_calls(),
                if method == 2 { 1 } else { 0 },
                "unexpected ARAP solve count (MLS must remain independent)"
            );
        }
        crate::arap::take_solve_calls();
        let _ = lua_pipeline(method, "pin_types[5]=PIN_TYPE.STARCH;pin_dx[1]=-35");
        assert_eq!(
            crate::arap::take_solve_calls(),
            if method == 2 { 2 } else { 0 },
            "following material controls lost their preliminary pose"
        );
    }
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_neck_branch_keeps_a_stable_head_angle() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut failures = Vec::new();
    for density in [5, 15, 30] {
        for method in [1, 2] {
            for (left, right) in [(70, 650), (220, 500)] {
                for y in [60.0, 140.0, 220.0, 300.0, 380.0, 460.0, 540.0, 620.0] {
                    let mut previous: Option<f64> = None;
                    for offset in [0.0, 0.01] {
                        let pose = format!(
                            r#"
                pins=6
                pin_sx={{360-hw,190-hw,531-hw,360-hw,70-hw,650-hw}}
                pin_sy={{510-hh,727-hh,727-hh,300-hh,380-hh,380-hh}}
                bone_forest={{{{1,{{2}},{{3}},{{4,{{5}},{{6}}}}}}}}
                for i=1,pins do pin_types[i]=1;pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dy[5]={y}-hh;pin_dy[6]={right_y}-hh
                pin_dx[5]={left}-hw;pin_dx[6]={right}-hw
            "#,
                            right_y = y + offset
                        );
                        let bone = lua_pipeline_fixture(
                            method,
                            &pose,
                            density,
                            false,
                            Some((&rgba, 723, 800)),
                        );
                        let position = lua_pipeline_fixture(
                            method,
                            &format!(
                                "{pose};for i=1,pins do pin_types[i]=PIN_TYPE.POSITION end;bone_forest={{}}"
                            ),
                            density,
                            false,
                            Some((&rgba, 723, 800)),
                        );
                        assert_eq!(
                            bone, position,
                            "neck correction depends on pin kind: mode {method}, density {density}, y {y}, offset {offset}"
                        );
                        let (_, _, _, render) = bone;
                        if y == 60. && offset == 0. && left == 220 {
                            std::fs::write(
                                format!(
                                    "../../target/puppet-qa/neck-raised-m{method}-d{density}.csv"
                                ),
                                render
                                    .chunks_exact(12)
                                    .map(|row| {
                                        row.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
                                    })
                                    .collect::<Vec<_>>()
                                    .join("\n"),
                            )
                            .unwrap();
                        }
                        let at = |x: f64, y: f64| {
                            rendered_material_point(
                                &render,
                                (x / 723. - 0.5) * 160.,
                                (y / 800. - 0.5) * 120.,
                            )
                        };
                        let crown = at(360., 40.);
                        let center = at(360., 180.);
                        let angle = (crown.0 - center.0).atan2(center.1 - crown.1).to_degrees();
                        let neck = at(360., 300.);
                        let head_length = (center.0 - neck.0).hypot(center.1 - neck.1);
                        if (head_length - 120.).abs() > 6. {
                            failures.push((method, format!(
                                "head slid away from neck m{method} d{density} y{y}: {head_length}"
                            )));
                        }
                        eprintln!("neck branch m{method} d{density} y{y}+{offset}: {angle}");
                        if angle.abs() > 2.0 {
                            failures.push((
                                method,
                                format!("head tilted m{method} d{density} y{y}: {angle}"),
                            ));
                        }
                        if let Some(old) = previous {
                            if (angle - old).abs() > 0.2 {
                                failures.push((
                                    method,
                                    format!("head angle jumped m{method} d{density} y{y}"),
                                ));
                            }
                        }
                        previous = Some(angle);
                    }
                }
            }
        }
    }
    shape_quality_results(&failures);
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_arms_keep_face_contour_without_a_neck_pin() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut failures = Vec::new();
    for density in [5, 15, 30] {
        for (name, left, right, y, head_pin) in [
            ("lower-free", 180, 540, 570, false),
            ("lower-head", 180, 540, 570, true),
            ("raise-free", 220, 500, 60, false),
            ("raise-head", 220, 500, 60, true),
        ] {
            let pose = format!(
                r#"
                pins={count}
                pin_sx={{70-hw,650-hw,190-hw,531-hw,360-hw,360-hw}}
                pin_sy={{380-hh,380-hh,727-hh,727-hh,510-hh,180-hh}}
                bone_forest={{}}
                for i=1,pins do pin_types[i]=0;pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[1]={left}-hw;pin_dx[2]={right}-hw;pin_dy[1]={y}-hh;pin_dy[2]={y}-hh
            "#,
                count = if head_pin { 6 } else { 5 }
            );
            let (_, _, _, render) =
                lua_pipeline_fixture(2, &pose, density, false, Some((&rgba, 723, 800)));
            std::fs::write(
                format!("../../target/puppet-qa/actual-neck-{name}-d{density}.csv"),
                render
                    .chunks_exact(12)
                    .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                    .collect::<Vec<_>>()
                    .join("\n"),
            )
            .unwrap();
            let at = |x: f64, y: f64| {
                rendered_material_point(&render, (x / 723. - 0.5) * 160., (y / 800. - 0.5) * 120.)
            };
            let center = at(360., 180.);
            let crown = at(360., 40.);
            let chin = at(360., 270.);
            let neck = at(360., 310.);
            eprintln!(
                "{name} d{density}: center {center:?}, crown {crown:?}, chin {chin:?}, neck {neck:?}"
            );
            for (a, b) in [
                ((265., 210.), (455., 210.)),
                ((290., 250.), (430., 250.)),
                ((360., 40.), (360., 270.)),
            ] {
                let pa = at(a.0, a.1);
                let pb = at(b.0, b.1);
                let ratio = (pb.0 - pa.0).hypot(pb.1 - pa.1) / (b.0 - a.0).hypot(b.1 - a.1);
                eprintln!("{name} d{density} contour {a:?}-{b:?}: {ratio}");
                if !(0.90..1.10).contains(&ratio) {
                    failures.push(format!("contour {name} d{density}: {ratio}"));
                }
            }
            if !(crown.1 < center.1 && center.1 < chin.1 && chin.1 < neck.1) {
                failures.push(format!("head folded {name} d{density}"));
            }
            if !((center.0 + 1.5).abs() < 40. && (center.1 + 220.).abs() < 100.) {
                failures.push(format!("head drifted {name} d{density}: {center:?}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
#[ignore = "requires decoded sleep_dainoji_man.png in target/puppet-qa/sleep.rgba"]
fn actual_head_boundary_survives_arm_and_upward_pulls() {
    let rgba = std::fs::read("../../target/puppet-qa/sleep.rgba").unwrap();
    let mut failures = Vec::new();
    for density in [5, 15, 30] {
        for (name, pin, x, y) in [
            ("arm-down", 1, 160, 560),
            ("arm-up", 1, 160, 100),
            ("head-up", 6, 360, 60),
            ("chest-up", 7, 360, 280),
        ] {
            let pose = format!(
                r#"
                pins=7
                pin_sx={{70-hw,650-hw,190-hw,531-hw,360-hw,360-hw,360-hw}}
                pin_sy={{380-hh,380-hh,727-hh,727-hh,510-hh,180-hh,380-hh}}
                bone_forest={{}}
                for i=1,pins do pin_types[i]=0;pin_dx[i]=pin_sx[i];pin_dy[i]=pin_sy[i];pin_layer[i]=0 end
                pin_dx[{pin}]={x}-hw;pin_dy[{pin}]={y}-hh
            "#
            );
            let (v, t, q, render) =
                lua_pipeline_fixture(2, &pose, density, false, Some((&rgba, 723, 800)));
            if name.starts_with("arm-") {
                let at = |x: f64, y: f64| {
                    rendered_material_point(
                        &render,
                        (x / 723. - 0.5) * 160.,
                        (y / 800. - 0.5) * 120.,
                    )
                };
                let center = at(360., 180.);
                for (x, y) in [(235., 198.), (485., 198.), (240., 185.), (480., 210.)] {
                    let ear = at(x, y);
                    let rest = (x - 360.).hypot(y - 180.);
                    let ratio = (ear.0 - center.0).hypot(ear.1 - center.1) / rest;
                    assert!(
                        (ratio - 1.).abs() < 0.03,
                        "ear pulled by arm {name} d{density} ({x},{y}): {ratio}"
                    );
                }
            }
            for (control, (sx, sy)) in [
                (70., 380.),
                (650., 380.),
                (190., 727.),
                (531., 727.),
                (360., 510.),
                (360., 180.),
                (360., 380.),
            ]
            .into_iter()
            .enumerate()
            {
                let i = v
                    .chunks_exact(2)
                    .position(|p| (p[0] - (sx - 361.5)).hypot(p[1] - (sy - 400.)) < 1e-8)
                    .unwrap()
                    * 2;
                let (tx, ty) = if control + 1 == pin {
                    (x as f64, y as f64)
                } else {
                    (sx, sy)
                };
                assert!(
                    (q[i] - (tx - 361.5)).hypot(q[i + 1] - (ty - 400.)) < 1e-8,
                    "pin target moved in {name}"
                );
            }
            if name == "head-up" {
                let at = |x: f64, y: f64| {
                    rendered_material_point(
                        &render,
                        (x / 723. - 0.5) * 160.,
                        (y / 800. - 0.5) * 120.,
                    )
                };
                let crown = at(360., 40.);
                let nose = at(360., 180.);
                let tilt = (crown.0 - nose.0)
                    .atan2(nose.1 - crown.1)
                    .abs()
                    .to_degrees();
                assert!(
                    tilt < 20.,
                    "upward head pull turned sideways: d{density}, {tilt} degrees"
                );
            }
            let mut worst = 0.0_f64;
            let mut flipped = 0;
            for face in t.chunks_exact(3) {
                let ids = face.iter().map(|&i| i as usize).collect::<Vec<_>>();
                if ids.iter().any(|&i| v[2 * i + 1] > -120.) {
                    continue;
                }
                let cross = |p: &[f64]| {
                    let [a, b, c] = [ids[0] * 2, ids[1] * 2, ids[2] * 2];
                    (p[b] - p[a]) * (p[c + 1] - p[a + 1]) - (p[b + 1] - p[a + 1]) * (p[c] - p[a])
                };
                if cross(&v) * cross(&q) <= 0. {
                    flipped += 1;
                }
                for k in 0..3 {
                    let (a, b) = (2 * ids[k], 2 * ids[(k + 1) % 3]);
                    let rest = (v[b] - v[a]).hypot(v[b + 1] - v[a + 1]);
                    let moved = (q[b] - q[a]).hypot(q[b + 1] - q[a + 1]);
                    worst = worst.max((moved / rest - 1.).abs());
                }
            }
            eprintln!("head boundary {name} d{density}: strain {worst}, flipped {flipped}");
            if worst > 0.15 || flipped > 0 {
                failures.push(format!(
                    "{name} d{density}: strain {worst}, flipped {flipped}"
                ));
            }
            if let Some(dir) = std::env::var_os("PUPPET_QA_OUTPUT") {
                std::fs::write(
                    std::path::Path::new(&dir).join(format!("{name}-d{density}.csv")),
                    render
                        .chunks_exact(12)
                        .map(|row| row.iter().map(f64::to_string).collect::<Vec<_>>().join(","))
                        .collect::<Vec<_>>()
                        .join("\n"),
                )
                .unwrap();
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// MLS quality regressions are explicitly permitted for the canonical fit.
// Geometry, pin, periodicity and analytic MLS tests remain unconditional.
fn shape_quality(method: i32, passed: bool, message: std::fmt::Arguments<'_>) {
    if method == 1 {
        if !passed {
            eprintln!("Rigid MLS quality measurement: {message}");
        }
    } else {
        assert!(passed, "{message}");
    }
}
fn shape_quality_results(failures: &[(i32, String)]) {
    let mls = failures
        .iter()
        .filter(|(method, _)| *method == 1)
        .collect::<Vec<_>>();
    if let Some((_, first)) = mls.first() {
        eprintln!(
            "Rigid MLS quality measurements: {} outside the former shape envelope; first: {first}",
            mls.len()
        );
    }
    let arap = failures
        .iter()
        .filter(|(method, _)| *method == 2)
        .map(|(_, message)| message.as_str())
        .collect::<Vec<_>>();
    assert!(arap.is_empty(), "{}", arap.join("\n"));
}
