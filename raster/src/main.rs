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


fn put_pixel(canvas: &render::Canvas<video::Window>, x: i32, y: i32, color: pixels::Color) {
	let canvas_x = x + HALF_CANVAS;
	let canvas_y = HALF_CANVAS - y - 1;
	canvas.pixel(canvas_x as i16, canvas_y as i16, color).unwrap();
}

fn interpolate(i0: i32, d0: i32, i1: i32, d1: i32) -> Vec<i32> {
	let mut values = vec![];
	let a: f32 = (d1 - d0) as f32 / (i1 - i0) as f32;
	let mut d = d0 as f32;
	for _ in i0..i1 {
		values.push(d as i32);
		d = d + a;
	}
	return values;
}

fn draw_filled_triangle(canvas: &render::Canvas<video::Window>, mut p0: rect::Point, mut p1: rect::Point, mut p2: rect::Point, color: pixels::Color) {
	// make p0 the lowest and p2 the highest point.
	if p1.y() < p0.y() { std::mem::swap(& mut p1, & mut p0); }
	if p2.y() < p0.y() { std::mem::swap(& mut p2, & mut p0); }
	if p2.y() < p1.y() { std::mem::swap(& mut p2, & mut p1); }

	// get horizontal points for each side
	let mut x01 = interpolate(p0.y(), p0.x(), p1.y(), p1.x());
	let mut x12 = interpolate(p1.y(), p1.x(), p2.y(), p2.x());
	let mut x02 = interpolate(p0.y(), p0.x(), p2.y(), p2.x());

	// two sides are x02 and (x01+x12)
	let x012 = [&x01[..], &x12[..]].concat();
	// remove duplicate y
	// x02.pop();

	// figure out which side is left
	let mut x_left = x012;
	let mut x_right = x02;
	let m = x_right.len() / 2;
	if x_right[m] < x_left[m] {
		std::mem::swap(& mut x_left, & mut x_right);
	}

	for y in p0.y()..p2.y() {
		for x in x_left[(y - p0.y()) as usize]..x_right[(y - p0.y()) as usize] {
			put_pixel(canvas, x, y, color);
		}
	}

}

fn draw_wireframe_triangle(canvas: &render::Canvas<video::Window>, mut p0: rect::Point, mut p1: rect::Point, mut p2: rect::Point, color: pixels::Color) {
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
		let ys = interpolate(p0.x(), p0.y(), p1.x(), p1.y());
		for x in p0.x()..p1.x() {
			put_pixel(canvas, x as i32, ys[(x - p0.x()) as usize] as i32, color);
		}
	} else {
		if p0.y() > p1.y() {
			let temp = p0;
			p0 = p1;
			p1 = temp;
		}
		let xs = interpolate(p0.y(), p0.x(), p1.y(), p1.x());
		for y in p0.y()..p1.y() {
			put_pixel(canvas, xs[(y-p0.y()) as usize], y as i32, color);
		}
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

	let p0 = rect::Point::new(-200, -250);
	let p1 = rect::Point::new(200, 50);
	let p2 = rect::Point::new(20, 250);
	let white = pixels::Color::RGB(255, 255, 255);
	let green = pixels::Color::RGB(0, 255, 0);

    // draw_line(&canvas, rect::Point::new(-50, -200), rect::Point::new(60, 240), pixels::Color::RGB(255, 255, 255));
    draw_filled_triangle(&canvas, p0, p1, p2, green);
    draw_wireframe_triangle(&canvas, p0, p1, p2, white);
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