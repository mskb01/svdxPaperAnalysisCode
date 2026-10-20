# svdxPaperAnalysisCode

Analysis code accompanying the SSANOVA-SVD paper.

## Organization

The `R/` directory contains analysis scripts for the nine datasets examined in the paper. Scripts prefixed with numbers `1`–`9` correspond to the respective datasets, with filenames indicating the dataset analyzed.

- `make*mainFigure.R` scripts generate the main figures for the corresponding datasets.
- `make*benchmark_supplement.R` scripts generate the supplementary benchmarking figures and results.

## Dependencies

The analysis code uses the [`svdxBenchmark`](https://github.com/mskb01/svdxBenchmark) repository.
