extern crate sdl2;

use sdl2::event::Event;
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Keycode;
use sdl2::pixels;
use sdl2::rect;
use sdl2::render;
use sdl2::video;

const SCREEN_WIDTH: i32 = 600;
const SCREEN_HEIGHT: i32 = 600;
const HALF_CANVAS: i32 = SCREEN_WIDTH / 2;
const VIEWPORT_WIDTH: f32 = 1.0;
const VIEWPORT_HEIGHT: f32 = 1.0;
const VIEWPORT_DEPTH: f32 = 1.0;

struct Pt {
    point: rect::Point,
    h: f32,
}

type Coords = [[f32; 4]; 4];

#[derive(Copy, Clone)]
struct Point3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Point3D {
    pub fn new(x: f32, y: f32, z: f32) -> Point3D {
        return Point3D { x: x, y: y, z: z };
    }
}

#[derive(Clone)]
struct Triangle {
    pub vertex1_idx: usize,
    pub vertex2_idx: usize,
    pub vertex3_idx: usize,
    pub color: pixels::Color,
}

impl Triangle {
    pub fn new(v1: usize, v2: usize, v3: usize, color: pixels::Color) -> Triangle {
        return Triangle {
            vertex1_idx: v1,
            vertex2_idx: v2,
            vertex3_idx: v3,
            color: color,
        };
    }
}

struct Model {
    vertexes: Vec<Point3D>,
    triangles: Vec<Triangle>,
    bounds_center: Point3D,
    bounds_radius: f32
}

struct Transform {
    scale: f32,
    rotation: Coords,
    translation: Point3D,
}

struct Instance<'a> {
    model: &'a Model,
    transform: Coords,
}

impl Instance<'_> {
    pub fn new(model: &Model, transform: Transform) -> Instance {
        let m = multiply_matrices(
            matrix_translate(transform.translation),
            multiply_matrices(transform.rotation, matrix_scale(transform.scale)),
        );
        // let m = multiply_matrices(
        // 	matrix_translate(transform.translation),
        // 	matrix_scale(1.0));
        // let m = multiply_matrices(matrix_translate(transform.translation), transform.rotation);
        return Instance {
            model: model,
            transform: m,
        };
    }
}

struct Camera {
    orientation: Coords,
    position: Point3D,
    clipping_planes: Vec<Plane>
}

struct Plane {
	normal: Point3D,
	distance: f32
}

fn put_pixel(canvas: &render::Canvas<video::Window>, x: i32, y: i32, color: pixels::Color) {
    let canvas_x = x + HALF_CANVAS;
    let canvas_y = HALF_CANVAS - y - 1;
    canvas
        .pixel(canvas_x as i16, canvas_y as i16, color)
        .unwrap();
}

fn interpolate(i0: i32, d0: f32, i1: i32, d1: f32) -> Vec<f32> {
    let mut values = vec![];
    let a: f32 = (d1 - d0) as f32 / (i1 - i0) as f32;
    let mut d = d0;
    for _ in i0..i1 {
        values.push(d);
        d = d + a;
    }
    return values;
}

fn matrix_rotation_y(degrees: f32) -> Coords {
    let radians = degrees * (std::f32::consts::FRAC_PI_2 / 180.0);
    let cos = radians.cos();
    let sin = radians.sin();
    return [
        [cos, 0.0, -sin, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [sin, 0.0, cos, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
}

fn matrix_scale(scale: f32) -> Coords {
    return [
        [scale, 0.0, 0.0, 0.0],
        [0.0, scale, 0.0, 0.0],
        [0.0, 0.0, scale, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
}

fn matrix_translate(t: Point3D) -> Coords {
    return [
        [1.0, 0.0, 0.0, t.x],
        [0.0, 1.0, 0.0, t.y],
        [0.0, 0.0, 1.0, t.z],
        [0.0, 0.0, 0.0, 1.0],
    ];
}

fn dot_product(v1: Point3D, v2: Point3D) -> f32 {
	v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

fn transpose_matrix(m: Coords) -> Coords {
    let mut result: Coords = [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] += m[j][i];
        }
    }
    return result;
}

fn multiply_matrices(m1: Coords, m2: Coords) -> Coords {
    let mut result: Coords = [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    for i in 0..4 {
        for j in 0..4 {
        	for k in 0..4 {
            	result[i][j] += m1[i][k] * m2[k][j];
        	}
        }
    }
    return result;
}

fn matrix_multiply_vector(v: &Point3D, m: Coords) -> Point3D {
    Point3D::new(
        v.x * m[0][0] + v.y * m[0][1] + v.z * m[0][2] + m[0][3],
        v.x * m[1][0] + v.y * m[1][1] + v.z * m[1][2] + m[1][3],
        v.x * m[2][0] + v.y * m[2][1] + v.z * m[2][2] + m[2][3],
    )
}

fn multiply_color(color: pixels::Color, factor: f32) -> pixels::Color {
    return pixels::Color::RGB(
        (color.r as f32 * factor) as u8,
        (color.g as f32 * factor) as u8,
        (color.b as f32 * factor) as u8,
    );
}

fn viewport_to_canvas(x: f32, y: f32) -> rect::Point {
    let canvas_x: i32 = (x * (SCREEN_WIDTH as f32 / VIEWPORT_WIDTH)) as i32;
    let canvas_y: i32 = (y * (SCREEN_HEIGHT as f32 / VIEWPORT_HEIGHT)) as i32;
    return rect::Point::new(canvas_x, canvas_y);
}

fn project_vertex(v: &Point3D) -> rect::Point {
    return viewport_to_canvas(v.x * VIEWPORT_DEPTH / v.z, v.y * VIEWPORT_DEPTH / v.z);
}

fn draw_wireframe_triangle(
    canvas: &render::Canvas<video::Window>,
    p0: rect::Point,
    p1: rect::Point,
    p2: rect::Point,
    color: pixels::Color,
) {
    draw_line(canvas, p0, p1, color);
    draw_line(canvas, p1, p2, color);
    draw_line(canvas, p2, p0, color);
}

fn draw_line(
    canvas: &render::Canvas<video::Window>,
    mut p0: rect::Point,
    mut p1: rect::Point,
    color: pixels::Color,
) -> () {
    if (p1.x() - p0.x()).abs() > (p1.y() - p0.y()).abs() {
        if p0.x() > p1.x() {
            let temp = p0;
            p0 = p1;
            p1 = temp;
        }
        let ys = interpolate(p0.x(), p0.y() as f32, p1.x(), p1.y() as f32);
        for x in p0.x()..p1.x() {
            put_pixel(canvas, x as i32, ys[(x - p0.x()) as usize] as i32, color);
        }
    } else {
        if p0.y() > p1.y() {
            let temp = p0;
            p0 = p1;
            p1 = temp;
        }
        let xs = interpolate(p0.y(), p0.x() as f32, p1.y(), p1.x() as f32);
        for y in p0.y()..p1.y() {
            put_pixel(canvas, xs[(y - p0.y()) as usize] as i32, y as i32, color);
        }
    }
}

fn render_triangle(
    canvas: &render::Canvas<video::Window>,
    triangle: &Triangle,
    projected: &Vec<rect::Point>,
) {
    draw_wireframe_triangle(
        canvas,
        projected[triangle.vertex1_idx],
        projected[triangle.vertex2_idx],
        projected[triangle.vertex3_idx],
        triangle.color,
    );
}

fn apply_transform(v: &Point3D, t: &Transform) -> Point3D {
	let v2 = v;
    // let v1 = Point3D::new(v.x * t.scale, v.y * t.scale, v.z * t.scale);
    // let v2 = matrix_multiply_vector(&v1, t.rotation);
    let v3 = Point3D::new(
        v2.x + t.translation.x,
        v2.y + t.translation.y,
        v2.z + t.translation.z,
    );
    return v3;
}

fn render_instance(canvas: &render::Canvas<video::Window>, model: &Model) -> () {
    let projected = model
        .vertexes
        .iter()
        .map(|v| {
            project_vertex(v)
        })
        .collect();

    for t in &model.triangles {
        render_triangle(canvas, t, &projected);
    }
}

fn clip_triangle(plane: &Plane, triangle: Triangle, vertexes: &Vec<Point3D>) -> Vec<Triangle> {
	let v1 = vertexes[triangle.vertex1_idx];
	let v2 = vertexes[triangle.vertex2_idx];
	let v3 = vertexes[triangle.vertex3_idx];

	// TODO: this seems different from the sphere clip test?
	let in1 = dot_product(plane.normal, v1) + plane.distance > 0.0;
	let in2 = dot_product(plane.normal, v2) + plane.distance > 0.0;
	let in3 = dot_product(plane.normal, v3) + plane.distance > 0.0;
	let in_count = in1 as i32 + in2 as i32 + in3 as i32;
	if in_count == 0 {
		// full clip, don't return anything.
		println!("clipped whole triangle");
		return vec![];
	} else if in_count == 3 {
		// preserve whole triangle
		println!("all points in plane");
		return vec![triangle];
	} else {
		println!("partial clipping possible, returning full triangle");
		return vec![triangle];
	}
}

fn transform_and_clip(planes: &Vec<Plane>, model: &Model, transform: Coords) -> Option<Model> {
	// early clip, only transform bounds center.

	let bounds_center = matrix_multiply_vector(&model.bounds_center, transform);
	let radius2 = model.bounds_radius * model.bounds_radius;
	for p in planes {
		// TODO: understand this. What does signed distance have to do with dot product?
		// get distance from center to plane.
		let distance2 = dot_product(p.normal, bounds_center) + p.distance;
		if distance2 < -radius2 {
			println!("Early discard !");
			return None;
		}
	}
    let positioned = model
        .vertexes
        .iter()
        .map(|v| {
			matrix_multiply_vector(v, transform)
        })
        .collect();

    let mut triangles = model.triangles.clone();
    for p in planes {
    	let mut new_triangles = vec![];
    	// clip the positioned triangles
    	for t in triangles {
    		new_triangles.append(&mut clip_triangle(p, t, &positioned));
    		// new_triangles.push(t);
    	}
    	triangles = new_triangles;
    }

    Some(Model {
    	triangles: triangles,
    	vertexes: positioned,
    	bounds_radius: model.bounds_radius,
    	bounds_center: model.bounds_center
    })
}

fn render_scene(
    canvas: &render::Canvas<video::Window>,
    camera: Camera,
    instances: Vec<Instance>,
) -> () {
    let rotation = transpose_matrix(camera.orientation);
    let pos = camera.position;
    let translation = matrix_translate(Point3D::new(-1.0 * pos.x, -1.0 * pos.y, -1.0 * pos.z));
    let camera_matrix = multiply_matrices(rotation, translation);

    for i in instances {
        let transform = multiply_matrices(camera_matrix, i.transform);
        /* camera.clipping_planes */
        match transform_and_clip(&camera.clipping_planes, i.model, transform) {
        	Some(m) => render_instance(canvas, &m),
        	None => {}
        }

    }
}

fn main() -> Result<(), String> {
    let sdl_context = sdl2::init()?;
    let video_subsys = sdl_context.video()?;
    let window = video_subsys
        .window(
            "computer graphics from scratch: raster",
            SCREEN_WIDTH as u32,
            SCREEN_HEIGHT as u32,
        )
        .position_centered()
        .opengl()
        .build()
        .map_err(|e| e.to_string())?;

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;

    canvas.set_draw_color(pixels::Color::RGB(0, 0, 0));
    canvas.clear();

    let vertexes = vec![
        Point3D::new(1.0, 1.0, 1.0),
        Point3D::new(-1.0, 1.0, 1.0),
        Point3D::new(-1.0, -1.0, 1.0),
        Point3D::new(1.0, -1.0, 1.0),
        Point3D::new(1.0, 1.0, -1.0),
        Point3D::new(-1.0, 1.0, -1.0),
        Point3D::new(-1.0, -1.0, -1.0),
        Point3D::new(1.0, -1.0, -1.0),
    ];

    // let white = pixels::Color::RGB(255, 255, 255);
    let red = pixels::Color::RGB(255, 0, 0);
    let green = pixels::Color::RGB(0, 255, 0);
    let blue = pixels::Color::RGB(0, 0, 255);
    let yellow = pixels::Color::RGB(255, 255, 0);
    let purple = pixels::Color::RGB(255, 0, 255);
    let cyan = pixels::Color::RGB(0, 255, 255);

    let triangles = vec![
        Triangle::new(0, 1, 2, red),
        Triangle::new(0, 2, 3, red),
        Triangle::new(4, 0, 3, green),
        Triangle::new(4, 3, 7, green),
        Triangle::new(5, 4, 7, blue),
        Triangle::new(5, 7, 6, blue),
        Triangle::new(1, 5, 6, yellow),
        Triangle::new(1, 6, 2, yellow),
        Triangle::new(4, 5, 1, purple),
        Triangle::new(4, 1, 0, purple),
        Triangle::new(2, 6, 7, cyan),
        Triangle::new(2, 7, 3, cyan),
    ];

    let cube = Model {
        vertexes: vertexes,
        triangles: triangles,

        // center should have bounding volume of sphere radius root 3.
        bounds_center: Point3D::new(0.0, 0.0, 0.0),
        bounds_radius: (3.0 as f32).sqrt() // (sqrt (1.0 + 2.0))
    };

    let cube1 = Instance::new(
        &cube,
        Transform {
            scale: 0.75,
            rotation: matrix_scale(1.0),
            translation: Point3D::new(-1.5, 0.0, 7.0),
        },
    );

    let cube2 = Instance::new(
        &cube,
        Transform {
            scale: 1.0,
            rotation: matrix_rotation_y(195.0),
            translation: Point3D::new(1.25, 2.5, 7.5),
        },
    );

    let s2: f32 = (2 as f32).sqrt();// 1.0 / (2 as f32).sqrt();
    let clipping_planes = vec![
    	Plane { normal: Point3D::new(0.0, 0.0, 0.0), distance: 1.0}, // near
    	Plane { normal: Point3D::new(s2, 0.0, s2), distance: 0.0}, // left,
    	Plane { normal: Point3D::new(-s2, 0.0, s2), distance: 0.0}, // right
    	Plane { normal: Point3D::new(0.0, -s2, s2), distance: 0.0}, // top
    	Plane { normal: Point3D::new(0.0, s2, s2), distance: 0.0}, // bottom
    ];

    let camera = Camera {
        position: Point3D::new(-3.0, 1.0, 2.0),
        orientation: matrix_rotation_y(-30.0),
        clipping_planes: clipping_planes
    };

    render_scene(&canvas, camera, vec![cube1, cube2]);

    canvas.present();

    println!("done");

    let mut events = sdl_context.event_pump()?;

    'main: loop {
        for event in events.poll_iter() {
            match event {
                Event::Quit { .. } => break 'main,

                Event::KeyDown {
                    keycode: Some(keycode),
                    ..
                } => {
                    if keycode == Keycode::Escape {
                        break 'main;
                    }
                }

                _ => {}
            }
        }
    }

    Ok(())
}
