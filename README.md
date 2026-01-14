# AA Endotyper - Local CNN Training

This workflow performs local CNN training for individual studies in the endotype analysis workflow.

## Overview

This workflow:
1. Takes a single study's data as input
2. Performs local CNN training on that study
3. Outputs trained models, visualizations, and processed data
4. Prepares outputs for downstream federated aggregation

## Usage

### Single Study Execution

```bash
nextflow run main.nf \
  --study_id SDY569 \
  --study_tidy_file data/SDY569_tidy.csv \
  --cpeptide_file data/SDY569_cpeptide_auc_tidy.csv
```

### Multiple Studies (Sequential Execution)

For processing multiple studies, you can run the workflow multiple times:

```bash
# Study 1
nextflow run main.nf \
  --study_id SDY569 \
  --study_tidy_file data/SDY569_tidy.csv \
  --cpeptide_file data/SDY569_cpeptide_auc_tidy.csv

# Study 2  
nextflow run main.nf \
  --study_id SDY797 \
  --study_tidy_file data/SDY797_tidy.csv \
  --cpeptide_file data/SDY797_cpeptide_auc_tidy.csv

# Study 3
nextflow run main.nf \
  --study_id SDY1737 \
  --study_tidy_file data/SDY1737_tidy.csv \
  --cpeptide_file data/SDY1737_cpeptide_auc_tidy.csv
```

## Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `study_id` | string | Yes | Unique identifier for the study (e.g., SDY569) |
| `study_tidy_file` | path | Yes | Path to the tidy feature CSV file |
| `cpeptide_file` | path | Yes | Path to the C-peptide AUC CSV file |

## Outputs

For each study, the workflow produces:

### Models
- `results/local_cnn/{study_id}/models/`: Trained CNN models in H5 format

### Results
- `results/local_cnn/{study_id}/*.csv`: 
  - `{study_id}_local_results.csv`: Combined results with predictions
  - `{study_id}_train.csv`: Training data for potential federated use
  - `{study_id}_test.csv`: Test data for potential federated use  
  - `{study_id}_metrics.csv`: Performance metrics

### Visualizations
- `results/local_cnn/{study_id}/figures/`: Training plots and prediction visualizations

## Integration with Federated Workflow

The outputs from this workflow are designed to be consumed by the federated aggregation workflow in `lifebit-platform-federated-aggregation`. The federated workflow will:

1. Collect outputs from multiple study runs
2. Perform federated CNN training using the prepared train/test splits
3. Aggregate results across all studies
