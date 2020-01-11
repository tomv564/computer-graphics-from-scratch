extern crate sdl2;

use sdl2::event::Event;
use sdl2::pixels;
use sdl2::rect;
use sdl2::keyboard::Keycode;
use sdl2::render;
use sdl2::video;
use sdl2::gfx::primitives::DrawRenderer;

const SCREEN_WIDTH: i32 = 600;
const SCREEN_HEIGHT: i32 = 600;
const HALF_CANVAS: i32 = SCREEN_WIDTH / 2;
const VIEWPORT_WIDTH: f32 = 1.0;
const VIEWPORT_HEIGHT: f32 = 1.0;
const VIEWPORT_DEPTH: f32 = 1.0;

struct Pt {
	point: rect::Point,
	h: f32
}

#[derive(Copy, Clone)]
struct Point3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Point3D {
    pub fn new(x: f32, y: f32, z: f32) -> Point3D {
    	return Point3D {x: x, y: y, z: z};
    }
}

struct Triangle {
	pub vertex1_idx: usize,
	pub vertex2_idx: usize,
	pub vertex3_idx: usize,
	pub color: pixels::Color
}

impl Triangle {
	pub fn new(v1: usize, v2: usize, v3: usize, color: pixels::Color) -> Triangle {
		return Triangle {
			vertex1_idx: v1,
			vertex2_idx: v2,
			vertex3_idx: v3,
			color: color
		}
	}
}

fn put_pixel(canvas: &render::Canvas<video::Window>, x: i32, y: i32, color: pixels::Color) {
	let canvas_x = x + HALF_CANVAS;
	let canvas_y = HALF_CANVAS - y - 1;
	canvas.pixel(canvas_x as i16, canvas_y as i16, color).unwrap();
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

fn multiply_color(color: pixels::Color, factor: f32) -> pixels::Color {
	return pixels::Color::RGB((color.r as f32 * factor) as u8, (color.g as f32 * factor) as u8, (color.b as f32 * factor) as u8);
}

fn viewport_to_canvas(x: f32, y: f32) -> rect::Point {
	let canvas_x: i32 = (x * (SCREEN_WIDTH as f32 / VIEWPORT_WIDTH)) as i32;
	let canvas_y: i32 = (y * (SCREEN_HEIGHT as f32 / VIEWPORT_HEIGHT)) as i32;
	return rect::Point::new(canvas_x, canvas_y);
}

fn project_vertex(v: &Point3D) -> rect::Point {
	return viewport_to_canvas(v.x * VIEWPORT_DEPTH / v.z, v.y * VIEWPORT_DEPTH / v.z);
}

fn draw_wireframe_triangle(canvas: &render::Canvas<video::Window>, p0: rect::Point, p1: rect::Point, p2: rect::Point, color: pixels::Color) {
	draw_line(canvas, p0, p1, color);
	draw_line(canvas, p1, p2, color);
	draw_line(canvas, p2, p0, color);
}

fn draw_line(canvas: &render::Canvas<video::Window>, mut p0: rect::Point , mut p1: rect::Point, color: pixels::Color) -> () {

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
			put_pixel(canvas, xs[(y-p0.y()) as usize] as i32, y as i32, color);
		}
	}
}

fn render_triangle(canvas: &render::Canvas<video::Window>, triangle: Triangle, projected: &Vec<rect::Point>) {
	draw_wireframe_triangle(canvas, projected[triangle.vertex1_idx], projected[triangle.vertex2_idx], projected[triangle.vertex3_idx], triangle.color);
}

fn render_object(canvas: &render::Canvas<video::Window>, vertexes: &Vec<Point3D>, triangles: Vec<Triangle>) -> () {
	let projected = vertexes.iter().map(|v| project_vertex(v)).collect();
	for t in triangles {
		render_triangle(canvas, t, &projected);
	}
}

fn main() -> Result<(), String> {
    let sdl_context = sdl2::init()?;
    let video_subsys = sdl_context.video()?;
    let window = video_subsys.window("computer graphics from scratch: raster", SCREEN_WIDTH as u32, SCREEN_HEIGHT as u32)
        .position_centered()
        .opengl()
        .build()
        .map_err(|e| e.to_string())?;

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;

    canvas.set_draw_color(pixels::Color::RGB(0, 0, 0));
    canvas.clear();

    let mut vertices = vec![
    	Point3D::new(1.0, 1.0, 1.0),
    	Point3D::new(-1.0, 1.0, 1.0),
    	Point3D::new(-1.0, -1.0, 1.0),
    	Point3D::new(1.0, -1.0, 1.0),
    	Point3D::new(1.0, 1.0, -1.0),
    	Point3D::new(-1.0, 1.0, -1.0),
    	Point3D::new(-1.0, -1.0, -1.0),
    	Point3D::new(1.0, -1.0, -1.0)
    ];

	let white = pixels::Color::RGB(255, 255, 255);
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

	for v in &mut vertices {
		v.x -= 1.5;
		v.z += 7.0;
	}

	render_object(&canvas, &vertices, triangles);

    canvas.present();

    let mut events = sdl_context.event_pump()?;

    'main: loop {
        for event in events.poll_iter() {

            match event {

                Event::Quit {..} => break 'main,

                Event::KeyDown {keycode: Some(keycode), ..} => {
                    if keycode == Keycode::Escape {
                        break 'main
                    }
                }

                _ => {}
            }
        }
    }

    Ok(())
}