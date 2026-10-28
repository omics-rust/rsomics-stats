# Frozen pasilla TMM oracle data

These are Git-only test assets. Both this directory and
`tests/normalization.rs` are excluded from the published Rust library package.
The Rust implementation's MIT OR Apache-2.0 license does not relabel these
data. No upstream statistical program implementation is included here.

## Data provenance and actual license context

`kept-counts.tsv` preserves the exact retained8066×7 count records from the
team's accepted original source workflow40606516204, not a fabricated or
selected mini-example. The workflow starts from14599×7 pasilla counts carried
by DESeq2 1.52.0 and uses the existing design-based filterByExpr defaults before
default weighted TMM. DESeq2 is a raw-data carrier, not a fitted analysis here.

The retained DESeq2 `DESCRIPTION` states `License: LGPL (>= 3)`; its exact
bytes are preserved as DESeq2-DESCRIPTION. The DESeq2 vignette states that
these extdata counts were copied from the pasilla package. No separate
per-data-file license declaration was located in the retained numeric files.
For this raw-derived count fixture, the upstream carrier's stated
LGPL-3.0-or-later context and notices are conservatively retained; this is not
a new MIT/Apache/CC0 license grant or a legal-clearance claim.

`library-state.tsv` is frozen numeric reference output from that source
workflow, using R4.6.1/Bioc3.23, edgeR4.10.7 and limma3.68.5. Running R or
edgeR does not, by itself, assign their program licenses to numeric output.
This notice records its distinct generated-output provenance without claiming
a newly established data license or rights clearance. See the
[GNU output FAQ](https://www.gnu.org/licenses/gpl-faq.en.html#GPLOutput).

The unmodified standard LGPL3/GPL3 texts were obtained from the authoritative
[LGPL text](https://www.gnu.org/licenses/lgpl-3.0.txt) and
[GPL text](https://www.gnu.org/licenses/gpl-3.0.txt). They retain that
conservative input-carrier licensing context, not an established license
for generated outputs or a relicensing of Rust code. LGPL3 incorporates
GPL3 terms, so both texts are retained.

Original-study citation: Brooks et al. (2011), *Conservation of an RNA
regulatory map between Drosophila and mammals*, Genome Research193–202,
[doi:10.1101/gr.108662.110](https://doi.org/10.1101/gr.108662.110), as recorded
in the retained DESeq2 bibliography. Data-carrier citation: Love, Huber and
Anders (2014), *Moderated estimation of fold change and dispersion for RNA-seq
data with DESeq2*, Genome Biology15:550,
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8),
from the carrier's CITATION metadata. Carrier authors/maintainer are preserved
in DESeq2-DESCRIPTION; [carrier repository](https://github.com/thelovelab/DESeq2).

## Frozen identity and interpretation

- DESeq2 1.52.0 source archive SHA256:
  `8c91699286336350e66eec132ce6fdf5bb4af78e2a4d015a5a61224f62a95984`.
- Original pasilla_gene_counts.tsv.gz SHA256:
  `4b5570008936786da365d110df1c83c8cffc646d6c667d6b15ecbfb3d9ee735d`.
- Source workflow witness SHA256:
  `091f476ca726487d23e83de7d62b5bffa7cea9021b533a4dda2043c89289a577`.
- Source manifest SHA256:
  `af682fed901644550411e6ae74e22b6bda45a7e567ca15e2dddc2c6b17027e97`.
- kept-counts.tsv:1298718bytes, SHA256
  `a200f55b099013e1b3d03a00ea20ee5fd6070699ff3d3f036885a1115dad92c8`.
- library-state.tsv:929bytes, SHA256
  `76984dcad0685f8c5334406a777dd8776012ebc6d89afe03e9c4e01b3e51ddd7`.
- DESeq2-DESCRIPTION SHA256:
  `686cbd9e39e77c7bd58889b2b87283836ff57ffeb8f3f5dbc45fb611d1e5867a`.

These exact UTF-8 TSV bytes, including quotes/spacing and original gene/sample
order, are copied without reformatting. No new source run, recalculation,
count subset or factor fitting generated this fixture.

The test supplies postfilter_library_size (column2 after the sample ID) and
checks TMM_factor (column3) under the original abs1e-10+rel1e-10 band. The
expected factors are from the original source output, not recomputed by the
function under test. Expected zero-based reference2 (untreated3) is the
accepted current native/source-algorithm regression value; the R source
library-state table does NOT emit a reference index and is not claimed to.

This one real-data case exercises non-unity normalization, reference choice,
ranking, weighted trimming and factor centering. A constant-unity shortcut or
reference0 cannot satisfy it. It does not establish full edgeR/limma workflow
compatibility, other normalization methods or biological conclusions.
