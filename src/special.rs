//! Distribution tails shared across the stats crates.
//!
//! `chi2_sf` is the chi-squared survival function P(X > x) for X ~ χ²(df),
//! matching `scipy.stats.chi2.sf`. Non-finite inputs are guarded before the
//! tail evaluation: `NaN` x yields `NaN`, `x <= 0` yields `1.0`, `x = +∞`
//! yields `0.0`, and a `df` that is `<= 0` or non-finite yields `NaN`.

use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Survival function of the chi-squared distribution, `P(X > x)`.
///
/// Value-exact to `scipy.stats.chi2.sf`; the survival function is evaluated
/// directly (not as `1 - cdf`) so it stays accurate in the far tail.
pub fn chi2_sf(x: f64, df: f64) -> f64 {
    if x.is_nan() || !df.is_finite() || df <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 1.0;
    }
    if x.is_infinite() {
        return 0.0;
    }
    match ChiSquared::new(df) {
        Ok(dist) => dist.sf(x),
        Err(_) => f64::NAN,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Golden values from scipy.stats.chi2.sf (scipy 1.17.1).
    const GOLDEN: &[(f64, f64, f64)] = &[
        (0.0, 1.0, 1.000_000_000_000_000_00e0),
        (0.5, 1.0, 4.795_001_221_869_533_70e-1),
        (1.0, 1.0, 3.173_105_078_629_111_51e-1),
        (3.333, 1.0, 6.790_291_340_185_901_36e-2),
        (3.841, 1.0, 5.001_368_376_395_669_28e-2),
        (6.635, 1.0, 9.999_419_574_042_536_44e-3),
        (10.0, 1.0, 1.565_402_258_002_548_17e-3),
        (0.5, 2.0, 7.788_007_830_714_048_78e-1),
        (2.0, 2.0, 3.678_794_411_714_424_45e-1),
        (5.991, 2.0, 5.001_161_502_657_908_85e-2),
        (1.0, 3.0, 8.012_519_569_012_008_79e-1),
        (7.815, 3.0, 4.999_390_297_488_388_17e-2),
        (11.34, 3.0, 1.002_251_761_691_246_20e-2),
        (0.1, 5.0, 9.998_376_833_880_774_36e-1),
        (11.07, 5.0, 5.000_961_862_240_545_92e-2),
        (20.0, 5.0, 1.249_730_563_031_377_27e-3),
        (50.0, 10.0, 2.669_083_424_904_495_13e-7),
        (15.0, 10.0, 1.320_618_562_877_205_47e-1),
        (2.5, 4.0, 6.446_357_929_354_278_32e-1),
        (30.0, 20.0, 6.985_366_069_940_986_13e-2),
        // deep tail, large df, and near-zero (sf≈1) to guard the boundaries.
        (200.0, 1.0, 2.088_487_583_762_568_84e-45),
        (120.0, 100.0, 8.440_668_109_369_188_49e-2),
        (1e-8, 1.0, 9.999_202_115_440_526_39e-1),
    ];

    #[test]
    fn matches_scipy() {
        for &(x, df, want) in GOLDEN {
            let got = chi2_sf(x, df);
            let tol = 1e-9 * want.abs().max(1e-300);
            assert!(
                (got - want).abs() <= tol,
                "chi2_sf({x}, {df}) = {got}, want {want}"
            );
        }
    }

    #[test]
    fn non_finite_guarded() {
        assert!(chi2_sf(f64::NAN, 1.0).is_nan());
        assert!(chi2_sf(1.0, f64::NAN).is_nan());
        assert_eq!(chi2_sf(-1.0, 1.0), 1.0);
        assert_eq!(chi2_sf(0.0, 1.0), 1.0);
        assert_eq!(chi2_sf(f64::INFINITY, 1.0), 0.0);
        assert!(chi2_sf(1.0, 0.0).is_nan());
        assert!(chi2_sf(1.0, -3.0).is_nan());
        // df = +inf must yield NaN, not panic (statrs gamma_ur panics on shape=inf).
        assert!(chi2_sf(1.0, f64::INFINITY).is_nan());
        assert!(chi2_sf(1.0, f64::NEG_INFINITY).is_nan());
    }
}
