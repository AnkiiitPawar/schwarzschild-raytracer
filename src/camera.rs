//! A static observer's pinhole camera.
//!
//! The camera's screen directions are measured in the observer's local
//! orthonormal frame. A local direction with radial part n_r and tangential
//! part n_⊥ corresponds to the coordinate velocity
//! ```text
//! v = n_r sqrt(1 - 2M/r) r̂ + n_⊥,
//! ```
//! because proper radial distance is dr / sqrt(1 - 2M/r). This choice also
//! gives the correct impact parameter b = r |n_⊥| / sqrt(1 - 2M/r).

use crate::geodesic::{Photon, MASS};
use crate::vec3::Vec3;

#[derive(Clone, Copy, Debug)]
pub struct Camera {
    pub pos: Vec3,
    pub forward: Vec3,
    pub right: Vec3,
    pub up: Vec3,
    pub tan_half_fov: f64,
    pub aspect: f64,
}

impl Camera {
    /// Camera at `distance`, looking at the black hole, at `inclination_deg`
    /// from the disk's axis (+y): 0° is face-on, 90° is edge-on.
    pub fn looking_at_origin(
        distance: f64,
        inclination_deg: f64,
        fov_deg: f64,
        aspect: f64,
    ) -> Self {
        let i = inclination_deg.to_radians();
        let pos = Vec3::new(0.0, i.cos(), -i.sin()) * distance;
        let forward = (-pos).normalized();
        let world_up = if i.sin().abs() < 1e-9 {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let right = forward.cross(world_up).normalized();
        let up = right.cross(forward);
        Self {
            pos,
            forward,
            right,
            up,
            tan_half_fov: (0.5 * fov_deg.to_radians()).tan(),
            aspect,
        }
    }

    /// Ray through screen coordinates (sx, sy) in [-1, 1]².
    /// Returns the photon and L_y/E of the *physical* photon (which travels
    /// from the scene into the camera, the opposite way to the traced ray).
    pub fn photon(&self, sx: f64, sy: f64) -> (Photon, f64) {
        let local = (self.forward
            + self.right * (sx * self.tan_half_fov * self.aspect)
            + self.up * (sy * self.tan_half_fov))
            .normalized();
        self.photon_with_local_direction(local)
    }

    pub fn photon_with_local_direction(&self, local: Vec3) -> (Photon, f64) {
        let r = self.pos.norm();
        let r_hat = self.pos * (1.0 / r);
        let n_r = local.dot(r_hat);
        let n_perp = local - r_hat * n_r;
        let lapse = (1.0 - 2.0 * MASS / r).sqrt();
        let vel = r_hat * (n_r * lapse) + n_perp;
        let l_over_e = self.pos.cross(vel) * (1.0 / lapse);
        (Photon { pos: self.pos, vel }, -l_over_e.y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodesic::{trace, Outcome, TraceSettings, B_CRITICAL};

    #[test]
    fn shadow_edge_matches_synges_formula() {
        // A static observer at r sees the shadow edge at sin α = (3√3 M / r) sqrt(1 - 2M/r).
        let r = 30.0;
        let cam = Camera::looking_at_origin(r, 90.0, 60.0, 1.0);
        let alpha_edge = (B_CRITICAL / r * (1.0 - 2.0 / r).sqrt()).asin();
        let settings = TraceSettings::default();
        for (factor, captured) in [(0.995, true), (1.005, false)] {
            let a = factor * alpha_edge;
            let dir = cam.forward * a.cos() + cam.right * a.sin();
            let (photon, _) = cam.photon_with_local_direction(dir);
            let hit_horizon = matches!(trace(photon, None, &settings), Outcome::Captured);
            assert_eq!(hit_horizon, captured, "angle factor {factor}");
        }
    }

    #[test]
    fn impact_parameter_of_a_camera_ray_matches_the_local_angle() {
        let r = 25.0;
        let cam = Camera::looking_at_origin(r, 60.0, 60.0, 1.0);
        let a: f64 = 0.2;
        let dir = cam.forward * a.cos() + cam.up * a.sin();
        let (photon, _) = cam.photon_with_local_direction(dir);
        let expected = r * a.sin() / (1.0 - 2.0 / r).sqrt();
        assert!((photon.impact_parameter() - expected).abs() < 1e-12);
    }
}
