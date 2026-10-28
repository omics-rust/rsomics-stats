use crate::{Result, StatsError};

/// Default TMM factors in caller sample order.
#[derive(Debug, Clone, PartialEq)]
pub struct TmmFactors {
    pub reference_sample: usize,
    pub factors: Vec<f64>,
}

/// Computes default weighted TMM on nonzero retained rows and supplied libraries.
///
/// Rows must be rectangular, finite and nonnegative; libraries must be positive,
/// finite, matching the sample width and have a finite total. Rows and samples
/// are never dropped or reordered. The reference index is zero-based and factors
/// have geometric mean one. See `docs/tmm.md` for the bounded upstream profile.
pub fn tmm_factors(retained_rows: &[&[f64]], libraries: &[f64]) -> Result<TmmFactors> {
    if retained_rows.is_empty() || libraries.is_empty() {
        return Err(StatsError::Empty);
    }
    let invalid = |reason| StatsError::InvalidTmmInput(reason);
    if libraries
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(invalid("library sizes must be finite and positive"));
    }
    if !libraries.iter().sum::<f64>().is_finite() {
        return Err(invalid("library size total must be finite"));
    }
    for row in retained_rows {
        if row.len() != libraries.len() {
            return Err(invalid("retained row width differs from library sizes"));
        }
        let mut nonzero = false;
        for &count in *row {
            if !count.is_finite() || count < 0.0 {
                return Err(invalid("retained counts must be finite and nonnegative"));
            }
            nonzero |= count > 0.0;
        }
        if !nonzero {
            return Err(invalid("retained rows must not be all zero"));
        }
    }
    let n = libraries.len();
    let quartiles: Vec<_> = (0..n)
        .map(|sample| {
            let values: Vec<_> = retained_rows.iter().map(|row| row[sample]).collect();
            quantile(&values, 0.75) / libraries[sample]
        })
        .collect();
    let reference = if quantile(&quartiles, 0.5) < 1e-20 {
        let masses: Vec<_> = (0..n)
            .map(|sample| {
                retained_rows
                    .iter()
                    .map(|row| row[sample].sqrt())
                    .sum::<f64>()
            })
            .collect();
        (1..n).fold(0, |best, sample| {
            if masses[sample] > masses[best] {
                sample
            } else {
                best
            }
        })
    } else {
        let mean = quartiles.iter().sum::<f64>() / n as f64;
        (1..n).fold(0, |best, sample| {
            if (quartiles[sample] - mean).abs() < (quartiles[best] - mean).abs() {
                sample
            } else {
                best
            }
        })
    };
    let reference_values: Vec<_> = retained_rows.iter().map(|row| row[reference]).collect();
    let mut observed = vec![0.0; retained_rows.len()];
    let mut factors = Vec::with_capacity(n);
    for sample in 0..n {
        for (value, row) in observed.iter_mut().zip(retained_rows) {
            *value = row[sample];
        }
        factors.push(tmm_pair(
            &observed,
            &reference_values,
            libraries[sample],
            libraries[reference],
        ));
    }
    let center = (factors.iter().map(|factor| factor.ln()).sum::<f64>() / n as f64).exp();
    for (sample, factor) in factors.iter_mut().enumerate() {
        *factor /= center;
        if !factor.is_finite() || *factor <= 0.0 {
            return Err(StatsError::InvalidTmmFactor {
                sample,
                value: *factor,
            });
        }
    }
    Ok(TmmFactors {
        reference_sample: reference,
        factors,
    })
}

fn tmm_pair(
    observed: &[f64],
    reference: &[f64],
    observed_library: f64,
    reference_library: f64,
) -> f64 {
    let mut ratios = Vec::with_capacity(observed.len());
    let mut abundance = Vec::with_capacity(observed.len());
    let mut variance = Vec::with_capacity(observed.len());
    for (&observed, &reference) in observed.iter().zip(reference) {
        let observed_fraction = observed / observed_library;
        let reference_fraction = reference / reference_library;
        let ratio = (observed_fraction / reference_fraction).log2();
        let absolute = (observed_fraction.log2() + reference_fraction.log2()) / 2.0;
        if ratio.is_finite() && absolute.is_finite() && absolute > -1e10 {
            ratios.push(ratio);
            abundance.push(absolute);
            variance.push(
                (observed_library - observed) / observed_library / observed
                    + (reference_library - reference) / reference_library / reference,
            );
        }
    }
    if ratios.is_empty() || ratios.iter().all(|ratio| ratio.abs() < 1e-6) {
        return 1.0;
    }
    let ratio_ranks = average_rank(&ratios);
    let abundance_ranks = average_rank(&abundance);
    let n = ratios.len() as f64;
    let lower_ratio = (n * 0.3).floor() + 1.0;
    let upper_ratio = n + 1.0 - lower_ratio;
    let lower_abundance = (n * 0.05).floor() + 1.0;
    let upper_abundance = n + 1.0 - lower_abundance;
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for index in 0..ratios.len() {
        if ratio_ranks[index] >= lower_ratio
            && ratio_ranks[index] <= upper_ratio
            && abundance_ranks[index] >= lower_abundance
            && abundance_ranks[index] <= upper_abundance
        {
            let value = ratios[index] / variance[index];
            let weight = 1.0 / variance[index];
            if !value.is_nan() {
                numerator += value;
            }
            if !weight.is_nan() {
                denominator += weight;
            }
        }
    }
    let mean = numerator / denominator;
    2.0_f64.powf(if mean.is_nan() { 0.0 } else { mean })
}

fn average_rank(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<_> = (0..values.len()).collect();
    order.sort_by(|&left, &right| values[left].total_cmp(&values[right]));
    let mut ranks = vec![0.0; values.len()];
    let mut first = 0;
    while first < order.len() {
        let mut end = first + 1;
        while end < order.len() && values[order[first]] == values[order[end]] {
            end += 1;
        }
        let rank = (first + 1 + end) as f64 / 2.0;
        for &index in &order[first..end] {
            ranks[index] = rank;
        }
        first = end;
    }
    ranks
}

fn quantile(values: &[f64], probability: f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let position = (sorted.len() - 1) as f64 * probability;
    let lower = position.floor() as usize;
    if lower + 1 == sorted.len() {
        sorted[lower]
    } else {
        sorted[lower] + (position - lower as f64) * (sorted[lower + 1] - sorted[lower])
    }
}
