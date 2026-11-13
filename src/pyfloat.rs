//! Python float semantics radiopyo relies on (`//`, `%`, `round(x, 6)`).

/// Python `a // b` for floats (CPython `float_divmod`).
pub fn floordiv(a: f64, b: f64) -> f64 {
    divmod(a, b).0
}

/// Python `a % b` for floats (CPython `float_divmod`).
pub fn pymod(a: f64, b: f64) -> f64 {
    divmod(a, b).1
}

/// CPython `float_divmod` (Objects/floatobject.c).
fn divmod(vx: f64, wx: f64) -> (f64, f64) {
    let mut m = vx % wx; // C fmod
    let mut div = (vx - m) / wx;
    if m != 0.0 {
        if (wx < 0.0) != (m < 0.0) {
            m += wx;
            div -= 1.0;
        }
    } else {
        m = 0.0f64.copysign(wx);
    }
    let floordiv = if div != 0.0 {
        let mut f = div.floor();
        if div - f > 0.5 {
            f += 1.0;
        }
        f
    } else {
        0.0f64.copysign(vx / wx)
    };
    (floordiv, m)
}

/// Python `round(x, 6)`.
pub fn round6(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    // Correctly rounded decimal formatting, as CPython's round(x, ndigits).
    format!("{x:.6}").parse().unwrap_or(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floordiv_matches_python() {
        assert_eq!(floordiv(7.0, 2.0), 3.0);
        assert_eq!(floordiv(-7.0, 2.0), -4.0);
        assert_eq!(floordiv(5.0, f64::INFINITY), 0.0);
        assert_eq!(floordiv(0.3, 0.1), 2.0);
    }

    #[test]
    fn pymod_matches_python() {
        assert_eq!(pymod(-7.0, 2.0), 1.0);
        assert_eq!(pymod(5.0, f64::INFINITY), 5.0);
        assert_eq!(pymod(7.0, 2.0), 1.0);
    }

    #[test]
    fn round6_matches_python() {
        assert_eq!(round6(0.1234565), 0.123456);
        assert_eq!(round6(2.0000004999), 2.0);
        assert_eq!(round6(0.5), 0.5);
    }
}
