use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::time::Instant;

use schwarzschild_raytracer::{render, RenderSettings};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // usage: schwarzschild-raytracer [width] [height] [inclination_deg] [output.png]
    let args: Vec<String> = std::env::args().collect();
    let defaults = RenderSettings::default();
    let settings = RenderSettings {
        width: args
            .get(1)
            .and_then(|v| v.parse().ok())
            .unwrap_or(defaults.width),
        height: args
            .get(2)
            .and_then(|v| v.parse().ok())
            .unwrap_or(defaults.height),
        inclination_deg: args
            .get(3)
            .and_then(|v| v.parse().ok())
            .unwrap_or(defaults.inclination_deg),
        ..defaults
    };
    let out = args.get(4).map(String::as_str).unwrap_or("blackhole.png");

    let start = Instant::now();
    let pixels = render(&settings);
    println!(
        "rendered {}x{} ({}x{} samples/pixel) in {:.1} s",
        settings.width,
        settings.height,
        settings.samples,
        settings.samples,
        start.elapsed().as_secs_f64()
    );

    if let Some(dir) = Path::new(out).parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir)?;
        }
    }
    let mut encoder = png::Encoder::new(
        BufWriter::new(File::create(out)?),
        settings.width as u32,
        settings.height as u32,
    );
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.write_header()?.write_image_data(&pixels)?;
    println!("wrote {out}");
    Ok(())
}
