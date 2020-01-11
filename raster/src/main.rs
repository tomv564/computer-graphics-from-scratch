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


struct Pt {
	point: rect::Point,
	h: f32
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

fn draw_shaded_triangle(canvas: &render::Canvas<video::Window>, pt0: Pt, pt1: Pt, pt2: Pt, color: pixels::Color) {
	let Pt { point: mut p0, h: h0 } = pt0;
	let Pt { point: mut p1, h: h1 } = pt1;
	let Pt { point: mut p2, h: h2 } = pt2;

	// make p0 the lowest and p2 the highest point.
	if p1.y() < p0.y() { std::mem::swap(& mut p1, & mut p0); }
	if p2.y() < p0.y() { std::mem::swap(& mut p2, & mut p0); }
	if p2.y() < p1.y() { std::mem::swap(& mut p2, & mut p1); }


	// get horizontal points for each side
	let x01 = interpolate(p0.y(), p0.x() as f32, p1.y(), p1.x() as f32);
	let h01 = interpolate(p0.y(), h0, p1.y(), h1);

	let x12 = interpolate(p1.y(), p1.x() as f32, p2.y(), p2.x() as f32);
	let h12 = interpolate(p1.y(), h1, p2.y(), h2);

	let x02 = interpolate(p0.y(), p0.x() as f32, p2.y(), p2.x() as f32);
	let h02 = interpolate(p0.y(), h0, p2.y(), h2);

	// two sides are x02 and (x01+x12)
	let x012 = [&x01[..], &x12[..]].concat();
	let h012 = [&h01[..], &h12[..]].concat();
	// remove duplicate y
	// x02.pop();

	// figure out which side is left
	let mut x_left = x012;
	let mut x_right = x02;
	let mut h_left = h012;
	let mut h_right = h02;
	let m = x_right.len() / 2;
	if x_right[m] < x_left[m] {
		std::mem::swap(& mut x_left, & mut x_right);
		std::mem::swap(& mut h_left, & mut h_right);
	}

	for y in p0.y()..p2.y() {
		let y_index = (y - p0.y()) as usize;
		let x_l = x_left[y_index] as i32;
		let x_r = x_right[y_index] as i32;
		let h_segment = interpolate(x_l, h_left[y_index], x_r, h_right[y_index]);
		for x in x_l .. x_r {
			let shaded_color = multiply_color(color, h_segment[(x - x_l) as usize]);
			put_pixel(canvas, x, y, shaded_color);
		}
	}

}

fn multiply_color(color: pixels::Color, factor: f32) -> pixels::Color {
	return pixels::Color::RGB((color.r as f32 * factor) as u8, (color.g as f32 * factor) as u8, (color.b as f32 * factor) as u8);
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
	let pt0 = Pt { point: p0, h: 0.3 };
	let pt1 = Pt { point: p1, h: 0.1 };
	let pt2 = Pt { point: p2, h: 1.0 };
	let white = pixels::Color::RGB(255, 255, 255);
	let green = pixels::Color::RGB(0, 255, 0);

    // draw_line(&canvas, rect::Point::new(-50, -200), rect::Point::new(60, 240), pixels::Color::RGB(255, 255, 255));
    draw_shaded_triangle(&canvas, pt0, pt1, pt2, green);
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