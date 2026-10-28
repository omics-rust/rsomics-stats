# Default TMM normalization

`normalization::tmm_factors` computes default weighted trimmed-mean-of-M-values
factors from borrowed retained count rows and caller-supplied library sizes.
Its reference index is zero-based; factors preserve the caller's sample order
and are centered to geometric mean one. It neither filters nor reorders rows,
recomputes libraries, nor constructs effective libraries or offsets.

Rows must be nonempty, rectangular, finite, nonnegative and not all zero.
Fractional counts are allowed. Libraries must match the row width, be positive
and finite, and have a finite total. Empty input, invalid input or nonpositive/
nonfinite centered factors return `StatsError`; no rows or samples are silently
removed. This intentionally excludes the upstream facade's all-zero-input,
implicit-library and configurable-method behavior.

The numerical profile follows edgeR4.10.7 default TMM: normalized75th-percentile
reference selection with its sparse-data fallback,30% M and5% A trimming,
average ranks for ties, weighted means and geometric centering. Only NaN terms
are omitted from weighted sums. There is no public rank/quantile API or
product-specific expression-analysis state.

## Consumers and provenance

Both `rsomics-edger/src/normalization.rs::prepare` and
`rsomics-limma/src/normalization.rs::prepare` need this exact profile. Their
filtering, named input policies and effective-library/offset states remain
local; edgeR QL and limma voom are not foundation responsibilities. A consumer
can collect borrowed retained-row slices without copying the count matrix;
that adapter allocation and the public input validation still have a cost.

The implementation reuses the team's identical current149-line kernels from
edgeR3a97e8f4ea510b6d83d9e3554b71e9eab0b37add and
limmaabb4a01cfc074caab53c47336953039f6ada2e60. Historical mathematical seed:
rsomics-tmm-norm3f37f84a8935c419a37378bc43dc09930da83067. Parser, matrix,
CLI and donor-specific finite-weight skipping were not retained.

Behavior authority is edgeR4.10.7 commit
27ba7b39d819e68df2ecdea9f32a4b6039b2a0b5, `R/normLibSizes.R` and
`R/calcNormFactors.R`. Those source bodies were read; upstream R implementation
text is not transplanted. This is a bounded default-TMM profile, not a claim
that all edgeR normalization methods or either complete analysis product are
qualified by the shared function.
