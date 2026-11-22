//! Black-hole ray tracer: integrates light rays backwards from a camera
//! through Schwarzschild spacetime and shades what they hit.

pub mod camera;
pub mod geodesic;
pub mod shading;
pub mod vec3;

use camera::Camera;
use geodesic::{trace, Disk, Outcome, TraceSettings};
use rayon::prelude::*;
use shading::{blackbody_rgb, flux_profile, redshift, star_field, tonemap};

#[derive(Clone, Copy, Debug)]
pub struct RenderSettings {
    pub width: usize,
    pub height: usize,
    /// Camera distance from the black hole, in units of M.
    pub distance: f64,
    /// Angle between the line of sight and the disk axis, degrees.
    pub inclination_deg: f64,
    pub fov_deg: f64,
    /// Samples per pixel along each axis (2 means 2×2 supersampling).
    pub samples: usize,
    pub disk: Disk,
    /// Emitted temperature at the hottest radius, kelvin.
    pub peak_temperature: f64,
    pub exposure: f64,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            width: 1280,
            height: 640,
            distance: 45.0,
            inclination_deg: 82.0,
            fov_deg: 48.0,
            samples: 2,
            disk: Disk {
                r_in: 6.0,
                r_out: 18.0,
            },
            peak_temperature: 3800.0,
            exposure: 1.4,
        }
    }
}

/// Linear-light colour seen along one camera ray.
pub fn shade(camera: &Camera, s: &RenderSettings, sx: f64, sy: f64) -> [f64; 3] {
    let (photon, b_axis) = camera.photon(sx, sy);
    let settings = TraceSettings {
        r_escape: 1.5 * s.distance,
        ..Default::default()
    };
    match trace(photon, Some(s.disk), &settings) {
        Outcome::Captured => [0.0; 3],
        Outcome::Escaped { direction } => star_field(direction),
        Outcome::Disk { hit, radius } => {
            let r_peak = s.disk.r_in * 49.0 / 36.0;
            let shape = flux_profile(radius, s.disk.r_in) / flux_profile(r_peak, s.disk.r_in);
            // Cosmetic only: faint trailing spiral banding so the disk reads as gas.
            let phi = hit.z.atan2(hit.x);
            let texture = 1.0 + 0.10 * (14.0 * radius.ln() + 2.0 * phi).sin();
            // Cosmetic only: fade the outer rim smoothly instead of a hard cut-off.
            let e = ((s.disk.r_out - radius) / (0.35 * s.disk.r_out)).clamp(0.0, 1.0);
            let fade = e * e * (3.0 - 2.0 * e);
            let g = redshift(radius, b_axis, s.distance);
            let t_emit = s.peak_temperature * shape.powf(0.25);
            // observed temperature scales with g, bolometric intensity with g⁴
            let colour = blackbody_rgb(g * t_emit);
            let intensity = g.powi(4) * shape * texture * fade;
            colour.map(|c| c * intensity)
        }
    }
}

/// Render to a row-major RGB8 buffer.
pub fn render(s: &RenderSettings) -> Vec<u8> {
    let aspect = s.width as f64 / s.height as f64;
    let camera = Camera::looking_at_origin(s.distance, s.inclination_deg, s.fov_deg, aspect);
    let mut pixels = vec![0u8; s.width * s.height * 3];
    pixels
        .par_chunks_mut(s.width * 3)
        .enumerate()
        .for_each(|(py, row)| {
            for px in 0..s.width {
                let mut acc = [0.0; 3];
                for a in 0..s.samples {
                    for b in 0..s.samples {
                        let fx = (px as f64 + (a as f64 + 0.5) / s.samples as f64) / s.width as f64;
                        let fy =
                            (py as f64 + (b as f64 + 0.5) / s.samples as f64) / s.height as f64;
                        let c = shade(&camera, s, 2.0 * fx - 1.0, 1.0 - 2.0 * fy);
                        for k in 0..3 {
                            acc[k] += c[k];
                        }
                    }
                }
                let n = (s.samples * s.samples) as f64;
                let rgb = tonemap(acc.map(|c| c / n), s.exposure);
                row[3 * px..3 * px + 3].copy_from_slice(&rgb);
            }
        });
    pixels
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_render_has_a_dark_shadow_and_a_bright_disk() {
        let s = RenderSettings {
            width: 64,
            height: 32,
            samples: 1,
            ..Default::default()
        };
        let img = render(&s);
        let pixel = |x: usize, y: usize| &img[3 * (y * s.width + x)..3 * (y * s.width + x) + 3];
        // just above centre: we look straight into the shadow
        assert!(pixel(32, 13).iter().all(|&c| c < 10));
        let brightest = img
            .chunks(3)
            .map(|p| p.iter().map(|&c| c as u32).sum::<u32>())
            .max()
            .unwrap();
        assert!(brightest > 300);
    }
}
