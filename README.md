# rsomics-stats

Library-only statistical components shared by rsomics products. Input and
analysis policies belong to the consuming product; this crate has no CLI.

- Multiple-testing adjustment and p-value combination.
- Welch, Mann–Whitney, Fisher and Hardy–Weinberg tests.
- Survival-function primitives.
- Default weighted TMM normalization over borrowed retained rows.

Fallible functions use `StatsError`; callers must follow each function's input
contract. Products retain their own identities, filtering, missing-value
mapping and analysis-state schemas.

```rust
use rsomics_stats::normalization::tmm_factors;

let rows: &[&[f64]] = &[&[10.0, 20.0], &[90.0, 180.0]];
let normalized = tmm_factors(rows, &[100.0, 200.0])?;
assert_eq!(normalized.factors, [1.0, 1.0]);
# Ok::<(), rsomics_stats::StatsError>(())
```

The bounded TMM profile, upstream provenance and concrete edgeR/limma
consumers are documented in [docs/tmm.md](docs/tmm.md). That shared numerical
contract does not make the two analysis products interchangeable or qualify
their complete statistical workflows.

Licensed under MIT OR Apache-2.0.

The normalization integration test and its oracle data under
`tests/fixtures/pasilla-tmm/` are Git-only and excluded from the published
crate. The fixtures retain scoped upstream data-carrier licensing notices;
see their `PROVENANCE.md`. They are not relabeled under the Rust code license.
