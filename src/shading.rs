//! Disk emission, relativistic redshift, colour and tone mapping.

use crate::geodesic::MASS;
use crate::vec3::Vec3;

/// Frequency ratio g = ν_observed / ν_emitted for light from a disk element
/// on a circular geodesic at radius `r` (angular velocity Ω = sqrt(M/r³)),
/// received by a static observer at `r_obs`. `b_axis` is the photon's
/// angular momentum about the disk axis divided by its energy.
///
/// ```text
/// g = sqrt(1 - 3M/r) / ( sqrt(1 - 2M/r_obs) (1 - Ω b_axis) )
/// ```
///
/// The numerator combines gravitational redshift and time dilation; the
/// denominator is the Doppler shift (blueshift on the approaching side).
pub fn redshift(r: f64, b_axis: f64, r_obs: f64) -> f64 {
    let omega = (MASS / (r * r * r)).sqrt();
    (1.0 - 3.0 * MASS / r).sqrt() / ((1.0 - 2.0 * MASS / r_obs).sqrt() * (1.0 - omega * b_axis))
}

/// Radial shape of a thin disk's emitted flux with a zero-torque inner edge at r_in:
/// F(r) ∝ (1 - sqrt(r_in / r)) / r³. The local temperature is T ∝ F^(1/4).
pub fn flux_profile(r: f64, r_in: f64) -> f64 {
    ((1.0 - (r_in / r).sqrt()) / (r * r * r)).max(0.0)
}

/// Approximate colour of a blackbody at `kelvin` as *linear* RGB (0..1),
/// valid for roughly 1,000-40,000 K. Uses a widely used piecewise fit to the
/// Planck locus (after T. Helland), whose output is display-encoded sRGB, so
/// it is decoded to linear light before any lighting maths.
pub fn blackbody_rgb(kelvin: f64) -> [f64; 3] {
    let t = kelvin.clamp(1000.0, 40_000.0) / 100.0;
    let r = if t <= 66.0 {
        255.0
    } else {
        329.698_727_446 * (t - 60.0).powf(-0.133_204_759_2)
    };
    let g = if t <= 66.0 {
        99.470_802_586_1 * t.ln() - 161.119_568_166_1
    } else {
        288.122_169_528_3 * (t - 60.0).powf(-0.075_514_849_2)
    };
    let b = if t >= 66.0 {
        255.0
    } else if t <= 19.0 {
        0.0
    } else {
        138.517_731_223_1 * (t - 10.0).ln() - 305.044_792_730_7
    };
    [r, g, b].map(|c| (c.clamp(0.0, 255.0) / 255.0).powf(2.2))
}

/// Sparse procedural star field on the celestial sphere.
pub fn star_field(dir: Vec3) -> [f64; 3] {
    let theta = dir.y.clamp(-1.0, 1.0).acos();
    let phi = dir.z.atan2(dir.x);
    let (cells_t, cells_p) = (2400.0, 4800.0);
    let i = (theta / std::f64::consts::PI * cells_t) as u64;
    let j = ((phi / std::f64::consts::TAU + 0.5) * cells_p) as u64;
    let mut h = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ j.wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    h ^= h >> 31;
    h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    h ^= h >> 29;
    let u = (h >> 11) as f64 / (1u64 << 53) as f64;
    if u < 0.0015 {
        let brightness = 0.15 + 1.5 * (u / 0.0015).powi(6);
        [brightness * 0.9, brightness * 0.95, brightness]
    } else {
        [0.0, 0.0, 0.0]
    }
}

/// Filmic exposure curve 1 - exp(-x * exposure) per channel, then gamma
/// encoding to 8 bits. Bright regions saturate towards white (hot gas),
/// dimmer regions keep their colour.
pub fn tonemap(c: [f64; 3], exposure: f64) -> [u8; 3] {
    c.map(|x| {
        let mapped = 1.0 - (-x.max(0.0) * exposure).exp();
        (mapped.powf(1.0 / 2.2) * 255.0).round() as u8
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn face_on_redshift_at_the_isco_is_sqrt_one_half() {
        assert!((redshift(6.0, 0.0, 1e15) - 0.5f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn approaching_side_is_blueshifted_relative_to_receding_side() {
        assert!(redshift(10.0, 8.0, 50.0) > redshift(10.0, -8.0, 50.0));
    }

    #[test]
    fn flux_vanishes_at_the_inner_edge_and_peaks_just_outside() {
        assert_eq!(flux_profile(6.0, 6.0), 0.0);
        let peak = (600..2000)
            .map(|k| k as f64 / 100.0)
            .max_by(|a, b| flux_profile(*a, 6.0).total_cmp(&flux_profile(*b, 6.0)))
            .unwrap();
        assert!((peak - 6.0 * 49.0 / 36.0).abs() < 0.02, "peak at {peak}"); // r = (49/36) r_in
    }

    #[test]
    fn hotter_blackbodies_are_bluer() {
        let warm = blackbody_rgb(3000.0);
        let hot = blackbody_rgb(15000.0);
        assert!(warm[2] / warm[0] < hot[2] / hot[0]);
    }
}
