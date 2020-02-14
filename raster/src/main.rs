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

// struct Pt {
//     point: rect::Point,
//     h: f32,
// }

type Coords = [[f32; 4]; 4];

type DepthBuffer = [f32; (SCREEN_HEIGHT*SCREEN_WIDTH) as usize];

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

    pub fn multiply(&self, f: f32) -> Point3D {
        return Point3D {x: self.x*f, y: self.y*f, z: self.z*f};
    }

    pub fn add(&self, v: &Point3D) -> Point3D {
        return Point3D {x: self.x + v.x, y: self.y + v.y, z: self.z + v.z};
    }
}

// impl Display for Point3D {
//     fn
// }

#[derive(Clone)]
struct Triangle {
    // pub vertex1_idx: usize,
    // pub vertex2_idx: usize,
    // pub vertex3_idx: usize,
    pub indexes: Vec<usize>,
    pub color: pixels::Color,
    pub normals: Vec<Point3D>
}

impl Triangle {
    pub fn new(v1: usize, v2: usize, v3: usize, color: pixels::Color, normals: Vec<Point3D>) -> Triangle {
        return Triangle {
            indexes: vec![v1, v2, v3],
            color: color,
            normals: normals
        };
    }
}

struct Model {
    vertexes: Vec<Point3D>,
    triangles: Vec<Triangle>,
    bounds_center: Point3D,
    bounds_radius: f32
}

fn generate_sphere(divs: i32, color: pixels::Color) -> Model {
    let mut vertexes: Vec<Point3D> = vec![];
    let mut triangles: Vec<Triangle> = vec![];

    let delta_angle = 2.0*std::f32::consts::FRAC_PI_2 / divs as f32;

    for d in 0..divs+1 {
        let y = (2.0 / divs as f32) * (d as f32 - (divs as f32)/2.0);
        let radius = (1.0 - y*y).sqrt();
        println!("generate sphere at height {}", y);
        for i in 0..divs {
            let x = radius * ((i as f32 * delta_angle).cos());
            let z = radius * ((i as f32 * delta_angle).sin());
            let vertex = Point3D::new(x, y, z);
            println!("adding vertex {} {} {}", vertex.x, vertex.y, vertex.z);
            vertexes.push(vertex);
        }
    }

    for d in 0..divs {
        for i in 0..divs - 1 {
            let i0 = d*divs + i;
            triangles.push(Triangle::new(i0 as usize, (i0+divs+1) as usize, (i0+1) as usize, color, vec![vertexes[i0 as usize], vertexes[(i0+divs+1) as usize], vertexes[(i0+1) as usize]]));
            triangles.push(Triangle::new(i0 as usize, (i0+divs) as usize, (i0+divs+1)  as usize, color, vec![vertexes[i0 as usize], vertexes[(i0+divs) as usize], vertexes[(i0+divs+1) as usize]]));
        }
    }

    Model {
        vertexes,
        triangles,
        bounds_center: Point3D::new(0.0, 0.0, 0.0),
        bounds_radius: 1.0
    }

}


struct Transform {
    scale: f32,
    rotation: Coords,
    translation: Point3D,
}

struct Instance<'a> {
    model: &'a Model,
    transform: Coords,
    orientation: Coords
}

impl Instance<'_> {
    pub fn new(model: &Model, transform: Transform) -> Instance {

        let m = multiply_matrices(
            &matrix_translate(&transform.translation),
            &multiply_matrices(&transform.rotation, &matrix_scale(transform.scale)),
        );
        // let m = multiply_matrices(
        // 	matrix_translate(transform.translation),
        // 	matrix_scale(1.0));
        // let m = multiply_matrices(matrix_translate(transform.translation), transform.rotation);
        return Instance {
            model: model,
            transform: m,
            orientation: transform.rotation
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

enum Light {
    Ambient { intensity: f32 },
    Point {intensity: f32, position: Point3D },
    Directional { intensity: f32, direction: Point3D}
}

// enum Shading {
//     Flat,
//     Gouraud,
//     Phong
// }

// struct Light {
//     light_type: LightType,
//     intensity: f32,
//     vector: ///Point3D
// }

// impl Light {
//     pub fn new(light_type: LightType, intensity: f32, vector: Point3D) -> Light {
//         Light { light_type, intensity, vector }
//     }
// }

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

fn matrix_translate(t: &Point3D) -> Coords {
    return [
        [1.0, 0.0, 0.0, t.x],
        [0.0, 1.0, 0.0, t.y],
        [0.0, 0.0, 1.0, t.z],
        [0.0, 0.0, 0.0, 1.0],
    ];
}

fn dot_product(v1: &Point3D, v2: &Point3D) -> f32 {
	v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

fn magnitude(v: &Point3D) -> f32 {
    dot_product(v, v).sqrt()
}

fn cross_product(v1: &Point3D, v2: &Point3D) -> Point3D {
    return Point3D::new(
        v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x
    )
}

fn transpose_matrix(m: &Coords) -> Coords {
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

fn multiply_matrices(m1: &Coords, m2: &Coords) -> Coords {
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

// fn draw_wireframe_triangle(
//     canvas: &render::Canvas<video::Window>,
//     p0: rect::Point,
//     p1: rect::Point,
//     p2: rect::Point,
//     color: pixels::Color,
// ) {
//     draw_line(canvas, p0, p1, color);
//     draw_line(canvas, p1, p2, color);
//     draw_line(canvas, p2, p0, color);
// }

fn edge_interpolate(y0: i32, v0: f32, y1: i32, v1: f32, y2: i32, v2: f32) -> (Vec<f32>, Vec<f32>) {
    // get horizontal points for each side
    let v01 = interpolate(y0, v0 as f32, y1, v1 as f32);
    let v12 = interpolate(y1, v1 as f32, y2, v2 as f32);
    let v02 = interpolate(y0, v0 as f32, y2, v2 as f32);

    // two sides are v02 and (v01+v12)
    let v012 = [&v01[..], &v12[..]].concat();

    return (v02, v012);
}

fn sorted_vertex_indices(indices: &Vec<usize>, projected: &Vec<rect::Point>) -> (usize, usize, usize) {
    let mut i0 = 0;
    let mut i1 = 1;
    let mut i2 = 2;

    // make p0 the lowest and p2 the highest point.
    if projected[indices[i1]].y() < projected[indices[i0]].y() { std::mem::swap(& mut i1, & mut i0); }
    if projected[indices[i2]].y() < projected[indices[i0]].y() { std::mem::swap(& mut i2, & mut i0); }
    if projected[indices[i2]].y() < projected[indices[i1]].y() { std::mem::swap(& mut i2, & mut i1); }

    (i0, i1, i2)
}

fn allowed_by_depth_buffer(db: &mut DepthBuffer, x: i32, y: i32, inv_z: f32) -> bool {
    let canvas_x = x + HALF_CANVAS;
    let canvas_y = HALF_CANVAS - y - 1;

    if canvas_x < 0 || canvas_x >= SCREEN_WIDTH || canvas_y < 0 || canvas_y >= SCREEN_HEIGHT {
        return false;
    }

    let offset = (canvas_x + SCREEN_WIDTH * canvas_y) as usize;
    if db[offset] < inv_z {
        if db[offset] > 0.0 {
            println!("overdraw at {}, {}, prev inv_z {} new inv_z {}", x, y, db[offset], inv_z);
        }
        (*db)[offset] = inv_z;
        return true;
    }
    return false;
}

fn compute_triangle_normal(v0: &Point3D, v1: &Point3D, v2: &Point3D) -> Point3D {
    // two lines on the same surface by subtracting v0 from v1 and v2.
    let v0v1: Point3D = v0.multiply(-1.0).add(v1);
    let v0v2: Point3D = v0.multiply(-1.0).add(v2);

    // perpendicular line to v0v1 and v0v2.
    return cross_product(&v0v1, &v0v2);
}

fn compute_lighting(point:Point3D, normal: Point3D, camera: &Camera, lights: &Vec<Light>) -> f32 {
    let mut illumination = 0.0;


    let directional_illumination = |light_vector: &Point3D, intensity: &f32| {        // diffuse lighting
        let mut lums = 0.0;

        // diffuse
        let cos_alpha = dot_product(light_vector, &normal) / (magnitude(light_vector) * magnitude(&normal));
        if cos_alpha > 0.0 {
            lums += cos_alpha * intensity
        }

        // specular
        let reflected = normal.multiply(2.0*dot_product(&normal, light_vector)).add(&light_vector.multiply(-1.0));
        let view = camera.position.add(&point.multiply(-1.0));
        let cos_beta = dot_product(&reflected, &view) / (magnitude(&reflected) * magnitude(&view));
        if cos_beta > 0.0 {
            let specular = 50;
            lums += cos_beta.powi(specular) * intensity;
        }

        return lums
    };

    for l in lights {
        match l {
            Light::Ambient {intensity} => {
                illumination += intensity;
            },
            Light::Directional {intensity, direction} => {
                let camera_matrix = transpose_matrix(&camera.orientation);
                let rotated_light = matrix_multiply_vector(direction, camera_matrix);
                illumination += directional_illumination(&rotated_light, intensity);
            }
            Light::Point {intensity, position } => {
                let camera_matrix = multiply_matrices(&transpose_matrix(&camera.orientation), &matrix_translate(&camera.position.multiply(-1.0)));
                let transformed_light = matrix_multiply_vector(position, camera_matrix);
                let light_vector = point.multiply(-1.0).add(&transformed_light);
                illumination += directional_illumination(&light_vector, intensity);
            }
        }
    }
    return illumination;
}

fn render_triangle(canvas: &render::Canvas<video::Window>, db: &mut DepthBuffer,
    triangle: &Triangle,
    vertexes: &Vec<Point3D>,
    projected: &Vec<rect::Point>,
    camera: &Camera,
    lights: &Vec<Light>,
    orientation: Coords
 ) {

    let (i0, i1, i2) = sorted_vertex_indices(&triangle.indexes, projected);
    let v0 = vertexes[triangle.indexes[i0]];
    let v1 = vertexes[triangle.indexes[i1]];
    let v2 = vertexes[triangle.indexes[i2]];

    // backface culling.
    // calculate triangle normal
    let normal = compute_triangle_normal(&vertexes[triangle.indexes[0]], &vertexes[triangle.indexes[1]], &vertexes[triangle.indexes[2]]);
    let to_triangle_center = vertexes[triangle.indexes[0]].add(&vertexes[triangle.indexes[0]]).add(&vertexes[triangle.indexes[0]]).multiply(-1.0/3.0);

    // - determine if angle from camera is greater than 90 degrees
    // if dot_product(&to_triangle_center, &normal) < 0.0 {
    //     println!("back face detected");
    //     return
    // }

    // flat shading: intensity from center
    // let center = Point3D::new(
    //     (v0.x + v1.x + v2.x)/3.0,
    //     (v0.y + v1.y + v2.y)/3.0,
    //     (v0.z + v1.z + v2.z)/3.0);
    // let intensity = compute_lighting(center, normal, camera, lights);


    let p0 = projected[triangle.indexes[i0]];
    let p1 = projected[triangle.indexes[i1]];
    let p2 = projected[triangle.indexes[i2]];

    let (x01, x012) = edge_interpolate(p0.y, p0.x as f32, p1.y, p1.x as f32, p2.y, p2.x as f32);
    // Note that we use the unprojected vertex Z-values here.
    let (iz01, iz012) = edge_interpolate(p0.y, 1.0 / v0.z, p1.y, 1.0 / v1.z, p2.y, 1.0 / v2.z);

    // use vertex normals
    let transform = multiply_matrices(&transpose_matrix(&camera.orientation), &orientation);
    let normal0 = matrix_multiply_vector(&triangle.normals[i0 as usize], transform);
    let normal1 = matrix_multiply_vector(&triangle.normals[i1 as usize], transform);
    let normal2 = matrix_multiply_vector(&triangle.normals[i2 as usize], transform);

    // gouraud shading: get lighting at triangle vertices:
    let l0 = compute_lighting(v0, normal0, camera, lights);
    let l1 = compute_lighting(v1, normal1, camera, lights);
    let l2 = compute_lighting(v2, normal2, camera, lights);

    println!("{} {} {} triangle vertices {} {} {} lighting {} {} {}", triangle.color.r, triangle.color.g, triangle.color.b, triangle.indexes[i0], triangle.indexes[i1], triangle.indexes[i2], l0, l1, l2);


    let (l01, l012) = edge_interpolate(p0.y, l0, p1.y, l1, p2.y, l2);


	// figure out which side is left
	let mut x_left = x01;
	let mut x_right = x012;
    let mut iz_left = iz01;
    let mut iz_right = iz012;
    let mut l_left = l01;
    let mut l_right = l012;
    println!("Interpolated edges lengths {} {}", x_left.len(), x_right.len());
    if (x_left.len() == 0) {
        println!("skipped!");
        return;
    }
	let m = x_right.len() / 2;
	if x_right[m] < x_left[m] {
		std::mem::swap(& mut x_left, & mut x_right);
        std::mem::swap(& mut iz_left, & mut iz_right);
        std::mem::swap(& mut l_left, & mut l_right);
	}

	for y in p0.y()..p2.y() {

        // for this line, get start and end of X and 1/Z
        let y_offset = (y  - p0.y()) as usize;
        let xl = x_left[y_offset];
        let xr = x_right[y_offset];
        let zl = iz_left[y_offset];
        let zr = iz_right[y_offset];
        let ll = l_left[y_offset];
        let lr = l_right[y_offset];

        // interpolate
        // println!("interpolate {} to {}", xl, xr);
        let z_scan = interpolate(xl as i32, zl, xr as i32, zr);
        let l_scan = interpolate(xl as i32, ll, xr as i32, lr);

		for x in xl as i32 ..xr as i32 {
            let x_index = (x - (xl as i32)) as usize;

            // println!("y {} x {} z_scan length {} index {}", y, x_index, z_scan.len(), x_index);
            let inv_z: f32 = z_scan[x_index];
            // println!("{} index {} inv_z {}", x, x_index, inv_z);

            if allowed_by_depth_buffer(db, x, y, inv_z) {
                put_pixel(canvas, x, y, multiply_color(triangle.color, l_scan[x_index as usize]));
            }
		}
	}
}

// fn apply_transform(v: &Point3D, t: &Transform) -> Point3D {
// 	let v2 = v;
//     // let v1 = Point3D::new(v.x * t.scale, v.y * t.scale, v.z * t.scale);
//     // let v2 = matrix_multiply_vector(&v1, t.rotation);
//     let v3 = Point3D::new(
//         v2.x + t.translation.x,
//         v2.y + t.translation.y,
//         v2.z + t.translation.z,
//     );
//     return v3;
// }

fn render_instance(canvas: &render::Canvas<video::Window>, depth_buffer: &mut DepthBuffer, model: &Model, camera: &Camera, lights: &Vec<Light>, orientation: Coords) -> () {
    let projected = model
        .vertexes
        .iter()
        .map(|v| {
            project_vertex(v)
        })
        .collect();

    for t in &model.triangles {
        render_triangle(canvas, depth_buffer, t, &model.vertexes, &projected, camera, lights, orientation);
    }
}

fn clip_triangle(plane: &Plane, triangle: Triangle, vertexes: &Vec<Point3D>) -> Vec<Triangle> {
	let v1 = vertexes[triangle.indexes[0]];
	let v2 = vertexes[triangle.indexes[1]];
	let v3 = vertexes[triangle.indexes[2]];

	// TODO: this seems different from the sphere clip test?
	let in1 = dot_product(&plane.normal, &v1) + plane.distance > 0.0;
	let in2 = dot_product(&plane.normal, &v2) + plane.distance > 0.0;
	let in3 = dot_product(&plane.normal, &v3) + plane.distance > 0.0;
	let in_count = in1 as i32 + in2 as i32 + in3 as i32;
	if in_count == 0 {
		// full clip, don't return anything.
		println!("clipped whole triangle");
		return vec![];
	} else if in_count == 3 {
		// preserve whole triangle
		// println!("all points in plane");
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
		let distance2 = dot_product(&p.normal, &bounds_center) + p.distance;
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
    camera: &Camera,
    instances: Vec<Instance>,
    lights: &Vec<Light>
) -> () {
    let rotation = transpose_matrix(&camera.orientation);
    let pos = camera.position;
    let translation = matrix_translate(&Point3D::new(-1.0 * pos.x, -1.0 * pos.y, -1.0 * pos.z));
    let camera_matrix = multiply_matrices(&rotation, &translation);
    let mut depth_buffer: DepthBuffer = [0.0; (SCREEN_HEIGHT*SCREEN_WIDTH) as usize];
    let orientation = instances[0 as usize].orientation;
    for i in instances {
        let transform = multiply_matrices(&camera_matrix, &i.transform);
        /* camera.clipping_planes */
        match transform_and_clip(&camera.clipping_planes, i.model, transform) {
        	Some(m) => render_instance(canvas, &mut depth_buffer, &m, camera, lights, orientation),
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
        Triangle::new(0, 1, 2, red,    vec![Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, 1.0)]),
        Triangle::new(0, 2, 3, red,    vec![Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, 1.0)]),
        Triangle::new(4, 0, 3, green,  vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 0.0, 0.0)]),
        Triangle::new(4, 3, 7, green,  vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 0.0, 0.0)]),
        Triangle::new(5, 4, 7, blue,   vec![Point3D::new(0.0, 0.0, -1.0), Point3D::new(0.0, 0.0, -1.0), Point3D::new(0.0, 0.0, -1.0)]),
        Triangle::new(5, 7, 6, blue,   vec![Point3D::new(0.0, 0.0, -1.0), Point3D::new(0.0, 0.0, -1.0), Point3D::new(0.0, 0.0, -1.0)]),
        Triangle::new(1, 5, 6, yellow, vec![Point3D::new(-1.0, 0.0, 0.0), Point3D::new(-1.0, 0.0, 0.0), Point3D::new(-1.0, 0.0, 0.0)]),
        Triangle::new(1, 6, 2, yellow, vec![Point3D::new(-1.0, 0.0, 0.0), Point3D::new(-1.0, 0.0, 0.0), Point3D::new(-1.0, 0.0, 0.0)]),
        Triangle::new(1, 0, 5, purple, vec![Point3D::new(0.0, 1.0, 0.0), Point3D::new(0.0, 1.0, 0.0), Point3D::new(0.0, 1.0, 0.0)]),
        Triangle::new(5, 0, 4, purple, vec![Point3D::new(0.0, 1.0, 0.0), Point3D::new(0.0, 1.0, 0.0), Point3D::new(0.0, 1.0, 0.0)]),
        Triangle::new(2, 6, 7, cyan,   vec![Point3D::new(0.0, -1.0, 0.0), Point3D::new(0.0, -1.0, 0.0), Point3D::new(0.0, -1.0, 0.0)]),
        Triangle::new(2, 7, 3, cyan,   vec![Point3D::new(0.0, -1.0, 0.0), Point3D::new(0.0, -1.0, 0.0), Point3D::new(0.0, -1.0, 0.0)])
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

    let sphere = generate_sphere(15, pixels::Color::RGB(0, 255, 0));
    let sphere1 = Instance::new(
        &sphere,
        Transform {
            scale: 1.5,
            rotation: matrix_scale(1.0),
            translation: Point3D::new(1.75, -0.5, 7.0)
        });


    let s2: f32 = (2 as f32).sqrt();// 1.0 / (2 as f32).sqrt();
    let clipping_planes = vec![
    	Plane { normal: Point3D::new(0.0, 0.0, 0.0), distance: 1.0}, // near
    	Plane { normal: Point3D::new(s2, 0.0, s2), distance: 0.0}, // left,
    	Plane { normal: Point3D::new(-s2, 0.0, s2), distance: 0.0}, // right
    	Plane { normal: Point3D::new(0.0, -s2, s2), distance: 0.0}, // top
    	Plane { normal: Point3D::new(0.0, s2, s2), distance: 0.0}, // bottom
    ];

    let camera = Camera {
        position: Point3D::new(-1.0, 1.0, -2.0),
        orientation: matrix_rotation_y(-30.0),
        clipping_planes: clipping_planes
    };

    let lights = vec![
        Light::Ambient { intensity: 0.2 },
        Light::Directional { intensity: 0.2, direction: Point3D::new(-1.0, 0.0, 1.0)},
        Light::Point { intensity: 0.2, position: Point3D::new(-3.0, 2.0, -10.0)}
    ];

    render_scene(&canvas, &camera, vec![cube1, cube2, sphere1], &lights);

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
