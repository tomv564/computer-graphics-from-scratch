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
	for i in i0..i1 {
		values.push(d as i32);
		d = d + a;
	}
	return values;
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
    draw_line(&canvas, rect::Point::new(-50, -200), rect::Point::new(60, 240), pixels::Color::RGB(255, 255, 255));
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