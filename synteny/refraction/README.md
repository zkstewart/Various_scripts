# refract_by_synteny.py

Projects ("refracts") GFF3 gene coordinates from one genome onto another genome's coordinate space using SyRI synteny predictions. This lets annotations from several genomes share one karyotype (e.g., a Circos density track), with syntenically equivalent regions approximately lined up.

The output is for **positional analyses only** (e.g., gene density). Feature lengths change when they are stretched or compressed, so do not use it to extract sequences.

## Usage

```bash
python3 refract_by_synteny.py -g input.gff3 -t target.fasta -s comparison.syri.out --input_side right -o refracted.gff3 --gaps_bed gaps.bed
```

| Option | Meaning |
|---|---|
| `-g` | GFF3 to refract |
| `-t` | Target genome lengths: a FASTA, a two-column TSV (`chrom<tab>length`) or a `.fai` |
| `-s` | SyRI `.out` file comparing the GFF3's genome with the target |
| `--input_side` | Which side of the SyRI file the GFF3's genome is on. SyRI puts the reference (`-r`) on the `left` and the query (`-q`) on the `right` |
| `-o` | Output GFF3 (default: stdout) |
| `-q` | Optional input genome lengths (default: inferred from the SyRI file) |
| `--translocations` | `ignore` (default), `intra` or `all`; see `trans` / `invtr` below |
| `--min_trans_length` | Minimum translocation size to use (default 10,000 bp) |
| `--keep_unmapped` | Keep features on chromosomes with no synteny to the target |
| `--gaps_bed` | Optional BED file of the interpolated gaps |
| `--min_gap_length` | Minimum gap size written to `--gaps_bed` (default 10,000 bp, input genome) |

## How positions are mapped

1. **Anchors:** SyRI's `SYNAL` (syntenic) and `INVAL` (inverted) alignments are the fixed points. Positions inside an anchor map linearly onto its target counterpart.
2. **Interpolation:** positions between anchors are spread linearly across the matching gap in the target. These gaps hold unaligned sequence, repeats, HDRs, translocations, duplications and indels. Chromosome ends act as anchors too, so every feature lands within the target chromosome.
3. **One decision per gene:** a gene and its child features (mRNA, exon, CDS and so on) share one mapping decision, made at the gene's midpoint, so a gene model is never split across modes.

If the warning *"SyRI reports coordinates ... exceeds its given length"* appears, `--input_side` is probably the wrong way round.

## Output GFF3: `refract_mode`

Every feature gets a `refract_mode=` attribute in column 9. Coordinates and chromosome names (column 1) are in the target genome, and `##sequence-region` lines describe the target.

| `refract_mode` | Meaning | Confidence |
|---|---|---|
| `syn` | The gene's midpoint lies in a syntenic alignment, so its position is anchored directly | High |
| `inv` | The gene's midpoint lies in an inversion. It is mapped in reverse and **its strand is flipped** | High |
| `interp` | The gene lies between anchors and its position is interpolated. Expect this for most genes, especially in repetitive clusters such as RGAs | Approximate; check the gap's `stretch_ratio` |
| `trans` | *Only with `--translocations intra/all`.* The gene lies in a translocated region and has been moved to the region's counterpart in the target, so it may appear out of order or on another chromosome (`all`) | Moderate |
| `invtr` | As `trans`, but the region is also inverted, so the strand is flipped | Moderate |
| `unmapped` | *Only with `--keep_unmapped`.* The gene's chromosome has no synteny to the target, so its coordinates are unchanged | None |

The script prints a count of genes (top-level features) per mode, and anchors per chromosome, to stderr.

## Gaps BED: `stretch_ratio`

Each line is one interpolated gap between anchors, given in **target** coordinates (BED: 0-based start, end-exclusive). A `#` header line comes first.

| Column | Content |
|---|---|
| 1–3 | Target chromosome, start, end |
| 4 | `input_region`: where the gap was in the input genome (`chrom:start-end`, 1-based) |
| 5 | `stretch_ratio` = target gap length ÷ input gap length |
| 6 | `log2_stretch_ratio`, useful as a colour scale centred on 0 |

Interpreting `stretch_ratio`:

| Value | Meaning | Effect on a density track |
|---|---|---|
| ≈ 1 | The gap is about the same size in both genomes | Little distortion |
| < 1 | **Compressed:** the input has more sequence here, e.g. an insertion or repeat expansion | `interp` genes are packed more tightly, so density is **inflated** |
| > 1 | **Stretched:** the target has more sequence here | `interp` genes are spread out, so density is **diluted** |
| 0 (`-inf`) | The gap shrinks to a point in the target (start = end) | All its genes collapse onto one position; Circos won't draw the gap unless you give it a minimum width |

A density peak made mostly of `interp` genes that sits in a gap with a low `stretch_ratio` (e.g. below 0.25, or `log2` below −2) is probably an artefact of compression, not a real cluster.

To use the BED as a Circos highlight track, drop the header line. If your karyotype uses different IDs, rename the chromosomes to match.
