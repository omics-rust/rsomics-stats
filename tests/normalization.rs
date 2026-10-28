use rsomics_stats::normalization::tmm_factors;

#[test]
fn retained_rows_produce_ordered_centered_factors_and_reject_invalid_domain() {
    let rows: &[&[f64]] = &[&[10.0, 20.0, 10.0, 20.0], &[90.0, 180.0, 90.0, 180.0]];
    let normalized = tmm_factors(rows, &[100.0, 200.0, 100.0, 200.0]).unwrap();
    assert_eq!(normalized.reference_sample, 0);
    assert_eq!(normalized.factors, [1.0; 4]);

    let number = |field: &str| field.trim_matches('"').trim().parse::<f64>().unwrap();
    let mut counts = include_str!("fixtures/pasilla-tmm/kept-counts.tsv").lines();
    let samples: Vec<_> = counts
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(|sample| sample.trim_matches('"'))
        .collect();
    let retained: Vec<Vec<_>> = counts
        .map(|row| row.split('\t').skip(1).map(number).collect())
        .collect();
    let source: Vec<Vec<_>> = include_str!("fixtures/pasilla-tmm/library-state.tsv")
        .lines()
        .skip(1)
        .map(|row| row.split('\t').collect())
        .collect();
    assert_eq!(retained.len(), 8066);
    assert_eq!(
        samples,
        source
            .iter()
            .map(|row| row[0].trim_matches('"'))
            .collect::<Vec<_>>()
    );
    let libraries: Vec<_> = source.iter().map(|row| number(row[2])).collect();
    let expected: Vec<_> = source.iter().map(|row| number(row[3])).collect();
    let borrowed: Vec<_> = retained.iter().map(Vec::as_slice).collect();
    let normalized = tmm_factors(&borrowed, &libraries).unwrap();
    assert_eq!(normalized.reference_sample, 2);
    assert_eq!(normalized.factors.len(), expected.len());
    for (actual, expected) in normalized.factors.iter().zip(expected) {
        assert!((actual - expected).abs() <= 1e-10 + 1e-10 * expected.abs());
    }

    for (rows, libraries) in [
        (&[][..], &[1.0][..]),
        (&[&[1.0, 2.0][..]][..], &[1.0][..]),
        (&[&[0.0, 0.0][..]][..], &[1.0, 1.0][..]),
        (&[&[1.0, f64::NAN][..]][..], &[1.0, 1.0][..]),
        (&[&[1.0, 2.0][..]][..], &[1.0, 0.0][..]),
    ] {
        assert!(tmm_factors(rows, libraries).is_err());
    }
}
