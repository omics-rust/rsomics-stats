/// Hardy-Weinberg equilibrium exact test (Wigginton et al. 2005, PMID:15789306).
///
/// Adapted from the algorithm in bcftools fill-tags.c (MIT licence), which
/// implements the Wigginton recursion. Treats each alt allele independently
/// against the reference, using only diploid biallelic counts.
///
/// # Arguments
/// * `nref` — number of reference allele copies among genotyped samples
///   (`2 * n_hom_ref + n_het`)
/// * `nalt` — number of alt allele copies (`2 * n_hom_alt + n_het`)
/// * `nhet` — count of heterozygous samples
///
/// # Returns
/// `(p_hwe, p_exc_het)` where:
/// * `p_hwe` — two-tailed HWE p-value (probability of observations as or more
///   unlikely than observed under HWE)
/// * `p_exc_het` — one-tailed excess-heterozygosity p-value (probability of
///   observing as many or more hets than seen)
///
/// Returns `(1.0, 1.0)` when there is no variation (`nref == 0 || nalt == 0`).
pub fn hwe_exact(nref: u32, nalt: u32, nhet: u32) -> (f64, f64) {
    if nref == 0 || nalt == 0 {
        return (1.0, 1.0);
    }
    let total = nref + nalt;
    let nrare = nref.min(nalt);

    // Midpoint: expected heterozygotes under HWE, adjusted to same parity as nrare.
    let mut mid = (f64::from(nrare) * f64::from(total - nrare) / f64::from(total)) as u32;
    if (nrare & 1) != (mid & 1) {
        mid += 1;
    }

    let n_probs = nrare as usize + 1;
    let mut probs = vec![0.0_f64; n_probs];

    // Start all probability mass at the midpoint, then recurse outward.
    probs[mid as usize] = 1.0;
    let mut sum = 1.0_f64;

    // Recurse downward from mid (fewer hets).
    let mut het = mid;
    let mut hom_r = (nrare - mid) / 2;
    let mut hom_c = (total / 2) - het - hom_r;
    while het >= 2 {
        let p = probs[het as usize] * f64::from(het) * f64::from(het - 1)
            / (4.0 * f64::from(hom_r + 1) * f64::from(hom_c + 1));
        probs[(het - 2) as usize] = p;
        sum += p;
        hom_r += 1;
        hom_c += 1;
        het -= 2;
    }

    // Recurse upward from mid (more hets).
    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = (total / 2) - het - hom_r;
    while het + 2 <= nrare {
        let p = probs[het as usize] * 4.0 * f64::from(hom_r) * f64::from(hom_c)
            / (f64::from(het + 2) * f64::from(het + 1));
        probs[(het + 2) as usize] = p;
        sum += p;
        hom_r = hom_r.saturating_sub(1);
        hom_c = hom_c.saturating_sub(1);
        het += 2;
    }

    // Normalize.
    for p in &mut probs {
        *p /= sum;
    }

    let obs = nhet as usize;
    let obs_p = if obs < n_probs { probs[obs] } else { 0.0 };

    // ExcHet: P(nhet >= observed) — one-tailed excess heterozygosity.
    let p_exc_het: f64 = probs[obs..].iter().sum();

    // HWE two-tailed: sum of all configurations as or less likely than observed.
    let p_hwe: f64 = probs
        .iter()
        .filter(|&&p| p <= obs_p + 1e-15)
        .sum::<f64>()
        .min(1.0);

    (p_hwe, p_exc_het)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn perfect_hwe_gives_high_p() {
        // 50 ref + 50 alt + 50 het → exactly HWE
        let (p_hwe, _) = hwe_exact(100, 100, 50);
        assert!(p_hwe > 0.9, "p_hwe={p_hwe}");
    }

    #[test]
    fn all_hom_alt_excess_hom() {
        // No hets at all — deviation from HWE
        let (p_hwe, p_exc_het) = hwe_exact(0, 200, 0);
        // nref==0 → trivial branch
        assert!(approx(p_hwe, 1.0, 1e-10));
        assert!(approx(p_exc_het, 1.0, 1e-10));
    }

    #[test]
    fn extreme_excess_het_gives_low_p_exc_het() {
        // 2 hom ref + 2 hom alt + 96 het (extreme excess)
        let nref = 2 * 2 + 96; // 100
        let nalt = 2 * 2 + 96; // 100
        let (_, p_exc_het) = hwe_exact(nref, nalt, 96);
        assert!(p_exc_het < 0.01, "p_exc_het={p_exc_het}");
    }

    #[test]
    fn no_variation_returns_one() {
        let (p, q) = hwe_exact(200, 0, 0);
        assert!(approx(p, 1.0, 1e-15));
        assert!(approx(q, 1.0, 1e-15));
    }
}
