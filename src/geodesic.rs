//! Null geodesics of the Schwarzschild metric (units G = c = M = 1, horizon at r = 2).
//!
//! Light rays stay in a plane, and in that plane the orbit u(φ) = 1/r obeys
//! ```text
//! d²u/dφ² + u = 3 M u².
//! ```
//! Treating (r, φ) as flat polar coordinates, this is exactly the path of a
//! particle under the central acceleration
//! ```text
//! x'' = -3 M h² x / r⁵,   h = |x × x'| (conserved),
//! ```
//! so rays can be integrated in 3D Cartesian coordinates with no coordinate
//! singularities and no need to choose the orbital plane by hand.

use crate::vec3::Vec3;

pub const MASS: f64 = 1.0;
pub const HORIZON: f64 = 2.0 * MASS;
pub const PHOTON_SPHERE: f64 = 3.0 * MASS;
/// Critical impact parameter: rays with b < 3√3 M are captured.
pub const B_CRITICAL: f64 = 5.196_152_422_706_632 * MASS;

#[derive(Clone, Copy, Debug)]
pub struct Photon {
    pub pos: Vec3,
    pub vel: Vec3,
}

impl Photon {
    /// Squared specific angular momentum h² = |x × v|², conserved along the ray.
    pub fn h2(&self) -> f64 {
        self.pos.cross(self.vel).norm2()
    }

    /// Impact parameter b = L/E. With the trajectory's conserved "energy"
    /// ε = v²/2 - M h²/r³ it is b = h / sqrt(2ε).
    pub fn impact_parameter(&self) -> f64 {
        let r = self.pos.norm();
        let h2 = self.h2();
        let eps = 0.5 * self.vel.norm2() - MASS * h2 / (r * r * r);
        (h2 / (2.0 * eps)).sqrt()
    }
}

fn accel(pos: Vec3, h2: f64) -> Vec3 {
    let r2 = pos.norm2();
    pos * (-3.0 * MASS * h2 / (r2 * r2 * r2.sqrt()))
}

/// One classical fourth-order Runge-Kutta step of size `dl`.
pub fn rk4_step(p: &Photon, h2: f64, dl: f64) -> Photon {
    let (x, v) = (p.pos, p.vel);
    let k1x = v;
    let k1v = accel(x, h2);
    let k2x = v + k1v * (0.5 * dl);
    let k2v = accel(x + k1x * (0.5 * dl), h2);
    let k3x = v + k2v * (0.5 * dl);
    let k3v = accel(x + k2x * (0.5 * dl), h2);
    let k4x = v + k3v * dl;
    let k4v = accel(x + k3x * dl, h2);
    Photon {
        pos: x + (k1x + k2x * 2.0 + k3x * 2.0 + k4x) * (dl / 6.0),
        vel: v + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (dl / 6.0),
    }
}

/// A geometrically thin disk in the plane y = 0.
#[derive(Clone, Copy, Debug)]
pub struct Disk {
    pub r_in: f64,
    pub r_out: f64,
}

#[derive(Clone, Copy, Debug)]
pub enum Outcome {
    /// Crossed the horizon (or circled near the photon sphere past the step limit).
    Captured,
    /// Hit the disk at this point.
    Disk { hit: Vec3, radius: f64 },
    /// Left the scene heading in this direction.
    Escaped { direction: Vec3 },
}

#[derive(Clone, Copy, Debug)]
pub struct TraceSettings {
    /// Step size as a fraction of the current radius.
    pub step_fraction: f64,
    /// Rays beyond this radius that are moving outwards have escaped.
    pub r_escape: f64,
    pub max_steps: usize,
}

impl Default for TraceSettings {
    fn default() -> Self {
        Self {
            step_fraction: 0.02,
            r_escape: 100.0,
            max_steps: 20_000,
        }
    }
}

pub fn trace(start: Photon, disk: Option<Disk>, settings: &TraceSettings) -> Outcome {
    let h2 = start.h2();
    let mut p = start;
    for _ in 0..settings.max_steps {
        let r = p.pos.norm();
        if r < HORIZON {
            return Outcome::Captured;
        }
        if r > settings.r_escape && p.pos.dot(p.vel) > 0.0 {
            return Outcome::Escaped {
                direction: p.vel.normalized(),
            };
        }
        let next = rk4_step(&p, h2, settings.step_fraction * r);
        if let Some(d) = disk {
            if p.pos.y * next.pos.y < 0.0 {
                // Crossed y = 0 during this step: interpolate the crossing point.
                let t = p.pos.y / (p.pos.y - next.pos.y);
                let hit = p.pos + (next.pos - p.pos) * t;
                let radius = hit.norm();
                if (d.r_in..=d.r_out).contains(&radius) {
                    return Outcome::Disk { hit, radius };
                }
            }
        }
        p = next;
    }
    Outcome::Captured
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_far(b: f64, x0: f64) -> Photon {
        Photon {
            pos: Vec3::new(-x0, b, 0.0),
            vel: Vec3::new(1.0, 0.0, 0.0),
        }
    }

    #[test]
    fn rays_inside_the_critical_impact_parameter_are_captured() {
        let s = TraceSettings {
            r_escape: 2000.0,
            ..Default::default()
        };
        assert!(matches!(
            trace(from_far(0.97 * B_CRITICAL, 1000.0), None, &s),
            Outcome::Captured
        ));
        assert!(matches!(
            trace(from_far(1.03 * B_CRITICAL, 1000.0), None, &s),
            Outcome::Escaped { .. }
        ));
    }

    #[test]
    fn weak_field_deflection_is_4m_over_b() {
        // Einstein's light-bending angle 4GM/(c² b), plus the next term 15πM²/(4b²).
        let b = 500.0;
        let x0 = 2.0e5;
        let s = TraceSettings {
            r_escape: 2.0 * x0,
            ..Default::default()
        };
        let Outcome::Escaped { direction } = trace(from_far(b, x0), None, &s) else {
            panic!("ray should escape");
        };
        let deflection = (-direction.y).atan2(direction.x);
        let expected = 4.0 * MASS / b + 15.0 * std::f64::consts::PI * MASS * MASS / (4.0 * b * b);
        // starting at finite distance misses about 2 * 2M/x0 of the total bending
        assert!(
            ((deflection - expected) / expected).abs() < 5e-3,
            "deflection {deflection}, expected {expected}"
        );
    }

    #[test]
    fn photon_sphere_is_a_circular_orbit() {
        let r = PHOTON_SPHERE;
        let mut p = Photon {
            pos: Vec3::new(r, 0.0, 0.0),
            vel: Vec3::new(0.0, 0.0, 1.0),
        };
        let h2 = p.h2();
        for _ in 0..2000 {
            p = rk4_step(&p, h2, 0.005);
        }
        // the orbit is unstable, but over ~half a revolution it stays put
        assert!((p.pos.norm() - r).abs() < 1e-6, "r = {}", p.pos.norm());
        assert!((p.impact_parameter() - B_CRITICAL).abs() < 1e-9);
    }
}
