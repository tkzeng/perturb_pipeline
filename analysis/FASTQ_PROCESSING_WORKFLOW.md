# FASTQ Processing Workflow

Complete workflow for `prepare_fastq.smk` with example data flows.

## Input Structure
```tsv
sample_id       fastq_dir
pool1:gex_1     /run1/fastq;/run2/fastq;/run3/fastq
pool1:guide_1   /run1/fastq;/run2/fastq;/run3/fastq
```

## Complete Processing Flow

### Step 1: `concatenate_lanes` 
**What it does:** Combines ALL raw FASTQ files for each sample across ALL runs and lanes

```
Sample: pool1:gex_1
fastq_dir: /run1/fastq;/run2/fastq;/run3/fastq

Input files found:
- /run1/fastq/gex_1_S1_L001_R1_001.fastq.gz
- /run1/fastq/gex_1_S1_L002_R1_001.fastq.gz  
- /run2/fastq/gex_1_S2_L001_R1_001.fastq.gz
- /run2/fastq/gex_1_S2_L002_R1_001.fastq.gz
- /run3/fastq/gex_1_S3_L001_R1_001.fastq.gz

Shell: cat all_files > concatenated_lanes/pool1:gex_1_R1.fastq.gz

Output: One concatenated file containing ALL sequencing data
Jobs: 1 per sample
```

### Step 2: `recover_barcodes_main_per_lane` (PARALLELIZED)
**What it does:** Processes each original raw lane file independently for barcode recovery

```
Parallel jobs for pool1:gex_1:
Job 1: /run1/fastq/gex_1_S1_L001_R1_001.fastq.gz → barcode_recovery/main/pool1:gex_1/run0/gex_1_S1_L001_recovered_R1.fastq.gz  
Job 2: /run1/fastq/gex_1_S1_L002_R1_001.fastq.gz → barcode_recovery/main/pool1:gex_1/run0/gex_1_S1_L002_recovered_R1.fastq.gz
Job 3: /run2/fastq/gex_1_S2_L001_R1_001.fastq.gz → barcode_recovery/main/pool1:gex_1/run1/gex_1_S2_L001_recovered_R1.fastq.gz
Job 4: /run2/fastq/gex_1_S2_L002_R1_001.fastq.gz → barcode_recovery/main/pool1:gex_1/run1/gex_1_S2_L002_recovered_R1.fastq.gz
Job 5: /run3/fastq/gex_1_S3_L001_R1_001.fastq.gz → barcode_recovery/main/pool1:gex_1/run2/gex_1_S3_L001_recovered_R1.fastq.gz

Process: Python barcode recovery script on individual lane files
Jobs: 5+ per sample (one per lane file across all runs)
```

### Step 3: `merge_barcode_recovery_main`
**What it does:** Combines all per-lane barcode recovery results for each sample

```
Input: All lane recovery files for pool1:gex_1
Shell: cat all_lane_recovered_files > barcode_recovery/main/pool1:gex_1_recovered_R1.fastq.gz

Output: Single recovered file per sample  
Jobs: 1 per sample
```

### Step 4: `recover_undetermined_per_lane` (PARALLELIZED)
**What it does:** Each job processes ONE undetermined lane file and SPLITS it by sample indices

```
Parallel jobs for pool1:
Job 1 processes: /run1/fastq/Undetermined_S0_L001_R1_001.fastq.gz
↓ Python script reads barcodes and splits by sample ↓
Outputs:
- undetermined_index/pool1/run0/Undetermined_S0_L001/AGGTCAGA+TATCTTGT_R1.fastq.gz (pool1:gex_1 reads)  
- undetermined_index/pool1/run0/Undetermined_S0_L001/CACTCAAT+CACAACTT_R1.fastq.gz (pool1:guide_1 reads)
- undetermined_index/pool1/run0/Undetermined_S0_L001/recovery_summary.txt

Job 2 processes: /run1/fastq/Undetermined_S0_L002_R1_001.fastq.gz  
↓ Same splitting process ↓
Similar per-sample outputs in separate directory...

Job 3 processes: /run2/fastq/Undetermined_S0_L001_R1_001.fastq.gz
↓ Same splitting process ↓  
Similar per-sample outputs in separate directory...

Process: Python script analyzes each lane's undetermined reads and splits them by sample barcode
Jobs: 5+ per pool (each job processes one undetermined lane file)
```

### Step 5: `create_undetermined_fastq` (Direct Lane Access)
**What it does:** Collects and concatenates all lane files for a specific sample

```
Input for pool1:gex_1: Dependencies on all lane recovery summaries
Process: 
  python create_undetermined_single_sample.py \
    --recovery-dirs undetermined_index/pool1/run0/Undetermined_S0_L001 \
                    undetermined_index/pool1/run0/Undetermined_S0_L002 \
                    undetermined_index/pool1/run1/Undetermined_S0_L001 \
                    undetermined_index/pool1/run1/Undetermined_S0_L002 \
                    undetermined_index/pool1/run2/Undetermined_S0_L001 \
    --sample-id pool1:gex_1

Script finds all files matching pool1:gex_1 barcode and concatenates them:
  AGGTCAGA+TATCTTGT_R1.fastq.gz from each lane directory

Output: undetermined_fastqs/pool1:gex_1_R1.fastq.gz (contains undetermined reads from ALL lanes)
Jobs: 1 per sample
```

### Step 6: `recover_barcodes_undetermined`
**What it does:** Additional barcode recovery on sample-specific undetermined reads

```
Input: undetermined_fastqs/pool1:gex_1_R1.fastq.gz
Process: Python barcode recovery script  
Output: barcode_recovery/undetermined/pool1:gex_1_recovered_R1.fastq.gz
Jobs: 1 per sample
```

### Step 7: `merge_all_fastqs`
**What it does:** Combines ALL processed file types for final merged output

```
Input for pool1:gex_1:
- concatenated_lanes/pool1:gex_1_R1.fastq.gz (main raw)
- barcode_recovery/main/pool1:gex_1_recovered_R1.fastq.gz (main recovered)  
- undetermined_fastqs/pool1:gex_1_R1.fastq.gz (undetermined raw)
- barcode_recovery/undetermined/pool1:gex_1_recovered_R1.fastq.gz (undetermined recovered)

Shell: cat all_four_files > merged_fastqs/pool1:gex_1_all_merged_R1.fastq.gz

Output: Final merged file containing ALL sequencing data for the sample
Jobs: 1 per sample
```

## File Types Available for Analysis

1. **main_raw**: `concatenated_lanes/{sample_id}_R1.fastq.gz` - Raw concatenated reads  
2. **main_recovered**: `barcode_recovery/main/{sample_id}_recovered_R1.fastq.gz` - Barcode-recovered main reads
3. **undetermined_raw**: `undetermined_fastqs/{sample_id}_R1.fastq.gz` - Sample-specific undetermined reads
4. **undetermined_recovered**: `barcode_recovery/undetermined/{sample_id}_recovered_R1.fastq.gz` - Barcode-recovered undetermined reads
5. **all_merged**: `merged_fastqs/{sample_id}_all_merged_R1.fastq.gz` - All four types combined