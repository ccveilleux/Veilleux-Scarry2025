# Warthog Genome Bioinformatics #

## 1. Download the SRA files (batched, based on list of SRRs): ##

#!/bin/bash
#SBATCH -N 2
#SBATCH -c 8
#SBATCH --mem=20G                   # total memory allocation
#SBATCH -t 0-15:00:00              # set a reasonable max runtime
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mem=16G        # Adjust memory if needed

module load sratoolkit-2.10.9-gcc-11.2.0

# Get the SRR ID for this array task
SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" QC_srr_list.txt)

# Download the SRA file
prefetch $SRR_ID --max-size 100G -O ./sra/

# Verify the file was downloaded
if [ -f ./sra/$SRR_ID/$SRR_ID.sra ]; then
    echo "$SRR_ID successfully downloaded."
else
    echo "Error: $SRR_ID failed to download."
    exit 1
fi

## 2. Convert SRA files to FASTQ format: ##
### performed in batches, changing list of SRRs each time ###

nano batch1_fasterq_dump.sh
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=40G                   # total memory allocation
#SBATCH -t 0-05:00:00              # set a reasonable max runtime
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mem=16G        # Adjust memory if needed

# Load SRA Toolkit
module load sratoolkit-2.10.9-gcc-11.2.0

# Paths
SRA_DIR=/scratch/cveille/warthogs/sra
FASTQ_DIR=/scratch/cveille/warthogs/fastq

# Make fastq output directory if needed
mkdir -p $FASTQ_DIR

# List of SRRs (Batch 1)
SRR_LIST=(
SRR19174515
SRR19174516
SRR19174517
SRR19174518
SRR19174519
SRR19174521
SRR19174522
SRR19174523
SRR19174524
SRR19174525
SRR19174526
SRR19174527
SRR19174529
SRR19174530
SRR19174532
)

for SRR in "${SRR_LIST[@]}"
do
  echo "Processing $SRR..."

  # Move into SRA subfolder
  cd ${SRA_DIR}/${SRR}

  # Run fasterq-dump
  fasterq-dump ${SRR}.sra --split-files --threads 8

  # Move fastq files to fastq directory
  mv ${SRR}_1.fastq ${SRR}_2.fastq $FASTQ_DIR/

  # Clean up to save space
  cd ..
  rm -r ${SRR}
done

## 3. Run FASTQC on the FASTQ files: ##

nano fastqc_batch.sh
--------
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8G                   # total memory allocation
#SBATCH -t 0-04:00:00              # set a reasonable max runtime
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mem=16G        # Adjust memory if needed
#SBATCH --array=0-9

#Load module
module load fastqc-0.12.1-gcc-11.2.0
# Set working directories
FASTQ_DIR=/scratch/cveille/warthogs/fastq
OUT_DIR=/scratch/cveille/warthogs/fastqc_results

# Make output directory if it doesn't exist
mkdir -p $OUT_DIR

# Run FastQC on all FASTQ files
fastqc --threads 4 -o $OUT_DIR $FASTQ_DIR/*.fastq


## 4. Run Trimmomatic on the FASTQ files ##
### performed in batches, changing list of SRRs each time ###

nano trimmomatic_batch.sh
--------
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=416G                   # total memory allocation
#SBATCH -t 0-02:00:00              # set a reasonable max runtime
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --array=0-9

module load trimmomatic-0.39-gcc-12.1.0

# Set directories
FASTQ_DIR=/scratch/cveille/warthogs/fastq/QCed
TRIM_DIR=/scratch/cveille/warthogs/fastq/trimmed
ADAPTERS=/scratch/cveille/warthogs/adapters.fa

mkdir -p $TRIM_DIR

# Define array of sample IDs
SAMPLES=(
SRR19174504
SRR19174505
SRR19174506
SRR19174507
SRR19174508
SRR19174510
SRR19174511
SRR19174512
SRR19174513
SRR19174514
)

# Select sample for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

trimmomatic PE -threads 4 \
  $FASTQ_DIR/${SAMPLE}_1.fastq $FASTQ_DIR/${SAMPLE}_2.fastq \
  $TRIM_DIR/${SAMPLE}_1_paired.fastq $TRIM_DIR/${SAMPLE}_1_unpaired.fastq \
  $TRIM_DIR/${SAMPLE}_2_paired.fastq $TRIM_DIR/${SAMPLE}_2_unpaired.fastq \
  ILLUMINACLIP:$ADAPTERS:2:30:10 \
  LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50


## 5. Run FastQC on the trimmed FASTQ files: ##
### using previous fastqc_batch.sh script ###


## 6. Download and index the reference genome: ##

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/906/955/GCF_016906955.1_ROS_Pafr_v1/GCF_016906955.1_ROS_Pafr_v1_genomic.fna.gz
gunzip GCF_016906955.1_ROS_Pafr_v1_genomic.fna.gz

nano index_ref_genome.sh
--------------
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G                   # total memory allocation
#SBATCH -t 0-5:00:00              # set a reasonable max runtime
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)

#Load modules
module load bwa-0.7.17-gcc-12.1.0
module load samtools-1.9-gcc-12.1.0

# Index for BWA
bwa index GCF_016906955.1_ROS_Pafr_v1_genomic.fna

# Index for samtools
samtools faidx GCF_016906955.1_ROS_Pafr_v1_genomic.fna

## 7. Align reads to the reference genome using BWA: ##

nano bwa_align.sh
---------
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G                   # total memory allocation
#SBATCH -t 0-8:00:00              # set a reasonable max runtime
#SBATCH --array=0-34
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --cpus-per-task=8

# Load necessary modules
module load bwa-0.7.17-gcc-12.1.0
module load samtools-1.9-gcc-12.1.0

# Define directories and reference genome
TRIM_DIR=/scratch/cveille/warthogs/fastq/trimmed
BAM_DIR=/scratch/cveille/warthogs/bam
GENOME=/scratch/cveille/warthogs/genome/GCF_016906955.1_ROS_Pafr_v1_genomic.fna

# Create BAM output directory if it doesn't exist
mkdir -p $BAM_DIR

# Define array of sample IDs
SAMPLES=(
SRR19174500
SRR19174504
SRR19174505
SRR19174506
SRR19174507
SRR19174508
SRR19174510
SRR19174511
SRR19174512
SRR19174513
SRR19174514
SRR19174515
SRR19174516
SRR19174517
SRR19174518
SRR19174519
SRR19174521
SRR19174522
SRR19174523
SRR19174524
SRR19174525
SRR19174526
SRR19174527
SRR19174529
SRR19174530
SRR19174532
SRR19174533
SRR19174534
SRR19174535
SRR19174537
SRR19174541
SRR19174542
SRR19174543
SRR19174551
SRR19174552
)

# Select sample based on the SLURM array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: $SAMPLE"

# Run BWA MEM to align paired-end reads
bwa mem -t 8 $GENOME $TRIM_DIR/${SAMPLE}_1_paired.fastq $TRIM_DIR/${SAMPLE}_2_paired.fastq > $BAM_DIR/${SAMPLE}.sam

# Convert SAM to BAM, sort, and index
samtools view -bS $BAM_DIR/${SAMPLE}.sam | samtools sort -o $BAM_DIR/${SAMPLE}_sorted.bam
samtools index $BAM_DIR/${SAMPLE}_sorted.bam

# Remove the intermediate SAM file to save space
rm $BAM_DIR/${SAMPLE}.sam

echo "Alignment and processing for sample $SAMPLE complete."

## 8. Mark and remove duplicates using Picard: ##

nano picard_mark_duplicates.sh
---------
#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G                   # total memory allocation
#SBATCH -t 0-4:00:00              # set a reasonable max runtime
#SBATCH --array=0-34
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=carrie.c.veilleux@gmail.com # send-to address

# Load necessary modules
module load picard-2.26.2-gcc-12.1.0

# Define directories
BAM_DIR=/scratch/cveille/warthogs/bam

# Define array of sample IDs
SAMPLES=(
SRR19174500
SRR19174504
SRR19174505
SRR19174506
SRR19174507
SRR19174508
SRR19174510
SRR19174511
SRR19174512
SRR19174513
SRR19174514
SRR19174515
SRR19174516
SRR19174517
SRR19174518
SRR19174519
SRR19174521
SRR19174522
SRR19174523
SRR19174524
SRR19174525
SRR19174526
SRR19174527
SRR19174529
SRR19174530
SRR19174532
SRR19174533
SRR19174534
SRR19174535
SRR19174537
SRR19174541
SRR19174542
SRR19174543
SRR19174551
SRR19174552
)


# Select sample based on SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: $SAMPLE"

# Input sorted BAM file (assumed to have been generated in the alignment step)
INPUT_BAM=${BAM_DIR}/${SAMPLE}_sorted.bam

# Output deduplicated BAM file and metrics file
OUTPUT_BAM=${BAM_DIR}/${SAMPLE}_dedup.bam
METRICS_FILE=${BAM_DIR}/${SAMPLE}_dedup.metrics.txt

# Run Picard MarkDuplicates
picard MarkDuplicates \
    INPUT=$INPUT_BAM \
    OUTPUT=$OUTPUT_BAM \
    METRICS_FILE=$METRICS_FILE \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

# Index the deduplicated BAM file
samtools index $OUTPUT_BAM

echo "Duplicate marking completed for sample $SAMPLE"

## 9. Check alignment statistics using samtools flagstat: ##
interactive -c 4 -N 1 -t 0-1:00

for bam in *_dedup.bam; do
    echo "Processing $bam"
    samtools flagstat "$bam" > "${bam%.bam}_flagstat.txt"
done

## 10. Generate a summary of the alignment statistics: ##
interactive -c 4 -N 1 -t 0-1:00

module load shpc/python/3.9.2-slim/module
nano extract_summary.py
---------
import glob
import csv
import re

# Prepare headers
header = ["Sample", "PercentDuplication", "Perc.Mapped", "Perc.ProperlyPaired", "Perc.Singletons", "Total.QC.Passed.Reads"]
data = {}

# Parse *_dedup.metrics.txt files
for filename in glob.glob("*_dedup.metrics.txt"):
    sample = filename.replace("_dedup.metrics.txt", "")
    with open(filename, "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("## METRICS CLASS") and i + 2 < len(lines):
                fields = lines[i + 2].strip().split("\t")  # This is the actual data line
                if len(fields) >= 9:
                    percent_dup = fields[8]
                    data[sample] = {"PercentDuplication": percent_dup}
                break

# Parse *_dedup_flagstat.txt files
for filename in glob.glob("*_dedup_flagstat.txt"):
    sample = filename.replace("_dedup_flagstat.txt", "")
    if sample not in data:
        data[sample] = {}

    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "in total" in line:
                total_qc_passed = int(line.split("+")[0].strip())
                data[sample]["Total.QC.Passed.Reads"] = total_qc_passed
            elif "mapped (" in line and "mate" not in line:
                match = re.search(r"\(([\d\.]+)%", line)
                if match:
                    data[sample]["Perc.Mapped"] = match.group(1)
            elif "properly paired" in line:
                match = re.search(r"\(([\d\.]+)%", line)
                if match:
                    data[sample]["Perc.ProperlyPaired"] = match.group(1)
            elif "singletons" in line:
                match = re.search(r"\(([\d\.]+)%", line)
                if match:
                    data[sample]["Perc.Singletons"] = match.group(1)

# Write output CSV
with open("combined_summary.csv", "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(header)
    for sample in sorted(data.keys()):
        row = [sample]
        for field in header[1:]:
            row.append(data[sample].get(field, "NA"))
        writer.writerow(row)

print("✅ Combined summary written to combined_summary.csv")

# Run the script to extract summary metrics
python extract_summary.py


## 11. Extract Exons 3, 4, and 5 of OPN1LW and index files: ##
interactive -c 4 -N 1 -t 0-1:00

#exon 3
for bamfile in *_dedup.bam; do
  samtools view -b "$bamfile" NC_062560.1:128939591-128939759 > "${bamfile%.bam}_exon3.bam"
done

#exon 4
for bamfile in *_dedup; do
  samtools view -b "$bamfile" NC_062560.1:128940730-128940892 > "${bamfile%.bam}_exon4.bam"
done

#exon 5
for bamfile in *_dedup; do
  samtools view -b "$bamfile" NC_062560.1:128942034-128942273 > "${bamfile%.bam}_exon5.bam"
done

for bam in *.bam; do samtools index "$bam"; done


## 12. Visualize the bam files in IGV ##