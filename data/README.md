# Data Directory

This directory is intended for processed count data files.

## Data Structure

When running analyses, the following data structure is expected:

```
data/
├── [dataset_name]/
│   ├── [dataset_name].txt      # Gene count matrix
│   └── ...
```

Each dataset count file should be a tab-delimited text file with:
- Column 1: GeneID (gene symbol or identifier)
- Remaining columns: Sample count data

## Sample Metadata

Sample metadata files are stored in the `config/` directory and describe:
- Sample identifiers
- Cell types or conditions
- Batch information
- Any other experimental factors

## Replicating Analyses

To replicate the analyses in this repository:

1. Download the raw or processed data from GEO accessions listed in Supplementary Table 1
2. Place processed count files in the appropriate subdirectories under `data/`
3. Ensure sample metadata in `config/` matches your data
4. Run the analysis scripts

## Reference Datasets

This project includes data from the following sources:
- Our laboratory samples (piPSC, hESC cultures)
- Published datasets (see Supplementary Table 1):
  - Choi et al. 2019
  - Choi et al. 2024
  - Yoshimatsu et al. 2021
  - Zhu et al. 2023
  - Xiang et al. 2024
  - Zhi et al. 2022
  - Selmi et al. 2021
  - Io et al. 2021
  - Secher et al. 2017

For questions about data access or structure, please open an issue on GitHub or contact the corresponding author.
