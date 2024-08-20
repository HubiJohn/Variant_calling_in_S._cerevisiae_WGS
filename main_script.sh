#!/bin/bash

#Tworzenie folderow
mkdir Projekt
mkdir -p Projekt/reads_raw
mkdir -p Projekt/qc_raw
mkdir -p Projekt/reads_trimmed
mkdir -p Projekt/qc_trimmed
mkdir -p Projekt/ref
mkdir -p Projekt/bams
mkdir -p Projekt/post_alignment
mkdir -p Projekt/bcftools
mkdir -p Projekt/freebayes


## Sciaganie danych

# funkcja walidujaca czy bioprojekt ma dobry identyfikator"
validate_project() {
	local project_id="$1"
	if [[ $project_id == PRJNA* ]]; then
    	echo "$project_id"
	else
    	echo "Invalid bioproject"
    	exit 1
	fi
}

# funkcja walidujaca czy run ma dobry identyfikator
# oraz sciagajaca podane runy
validate_sra() {
	local sra_ids=("$@")
	local valid=true
	for sra_id in "${sra_ids[@]}"; do
    	if [[ $sra_id == SRR* ]]; then
        	echo "$sra_id"
        	fastq-dump -X 1000000 "$sra_id" \
                --split-files -O ./Projekt/reads_raw
    	else
        	echo "Invalid SRA ID: $sra_id"
        	valid=false
    	fi
	done
	if [ "$valid" = false ]; then
    	exit 1
	fi
}

# Rozprawienie sie z podanymi argumentami - podany projtekt czy runy?
while [[ "$#" -gt 0 ]]; do
	case "$1" in
    	-p|--project)
        	if [ -n "$2" ] && [[ "$2" != -* ]]; then
            	project_id="$2"
            	shift 2
        	else
            	echo "Error: Argument for $1 is missing" >&2
            	exit 1
        	fi
        	;;
    	-r|--sra)
        	while [ -n "$2" ] && [[ "$2" != -* ]]; do
            	sra_ids+=("$2")
            	shift
        	done
        	shift
        	;;
    	*)
        	echo "Invalid option: $1" >&2
        	exit 1
        	;;
	esac
done


# Sprawdzamy czy podany jest run czy projrkt i walidujemy
if [ -n "$project_id" ]; then
	validate_project "$project_id"   
fi


if [ -n "${sra_ids+set}" ]; then
	validate_sra "${sra_ids[@]}"
fi


# jezeli podany jest projekt, sciagnij wszystkie runy z niego 
if [ -n "$project_id" ]; then
	esearch -db sra -query "$project_number" \
    	| efetch -format runinfo -mode xml \
    	| xtract -pattern SraRunInfo -element Run \
    	| awk '{for (i=1; i<=NF; i++) print $i}' > SRR_list.txt
   	 
	while IFS= read -r run_number; do
    	fastq-dump -X 1000000 "$run_number" \
            --split-files -O ./Projekt/reads_raw
	done < SRR_list.txt
fi


## analiza jakosci 1

for raw_read in ./Projekt/reads_raw/*.fastq; do
	fastqc -o ./Projekt/qc_raw "$raw_read"
done

multiqc ./Projekt/qc_raw -o ./Projekt/qc_raw/







## edycja danych (trimming)

cat > Projekt/ref/TrueSeq_adapters.fasta <<EOF
>Read1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
EOF

for raw_read in Projekt/reads_raw/*_1.fastq; do
    run_id=$(basename "$raw_read" _1.fastq)
    trimmomatic PE \
		"./Projekt/reads_raw/${run_id}_1.fastq"\
            "./Projekt/reads_raw/${run_id}_2.fastq" \
		"./Projekt/reads_trimmed/${run_id}_1.fastq" \
            "./Projekt/reads_trimmed/${run_id}_1_unpaired_1.fa" \
		"./Projekt/reads_trimmed/${run_id}_2.fastq" \
            "./Projekt/reads_trimmed/${run_id}_2_unpaired.fa" \
		SLIDINGWINDOW:4:30 \
        ILLUMINACLIP:Projekt/ref/TruSeq_adapters.fasta:2:30:10 \
        MINLEN:70
done



## analiza jakości 2 (po edycji)


for raw_read in ./Projekt/reads_trimmed/*.fastq; do
	fastqc -o ./Projekt/qc_trimmed "$raw_read"
done

multiqc ./Projekt/qc_trimmed -o ./Projekt/qc_trimmed/




## mapowanie

wget -O Projekt/ref/R64-1-1.fasta https://ftp.ensembl.org/pub/release-112/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

bwa index Projekt/ref/R64-1-1.fasta

for read1 in Projekt/reads_trimmed/*1.fastq; do
    read2="${read1/1.fastq/2.fastq}"
    run_id=$(basename "$read1" "_1.fastq")
    OUTPUT_PREFIX="${READS_DIR}/${run_id}aligned"
    bwa mem Projekt/ref/R64-1-1.fasta $read1 $read2 \
        | samtools view -b -o "Projekt/bams/${run_id}.bam"
done




## analiza po mapowaniu (post-alignment)

for bam_file in Projekt/bams/*.bam; do
  run_id=$(basename "$bam_file" .bam)
  flagstat_output="Projekt/post_alignment/${run_id}_flagstat.txt"
  depth_output="Projekt/post_alignment/${run_id}_depth.txt"

  # indexowanie
  samtools sort -o "${bam_file}.tmp" "$bam_file" \
    && mv "${bam_file}.tmp" "$bam_file"
  samtools index "$bam_file"

  # fixmate
  samtools sort -n -o "${bam_file}.querysorted" "$bam_file" \
    && mv "${bam_file}.querysorted" "$bam_file"
  samtools fixmate -m "$bam_file" "${bam_file}.tmp" \
    && mv "${bam_file}.tmp" "$bam_file"

  # markdup
  samtools sort -o "${bam_file}.tmp" "$bam_file" \
    && mv "${bam_file}.tmp" "$bam_file"
  samtools markdup "$bam_file" "${bam_file}.tmp" \
    && mv "${bam_file}.tmp" "$bam_file"

  # flagstats i depth
  samtools flagstat "$bam_file" > "$flagstat_output"
  samtools depth "$bam_file" > "$depth_output"
done




## variant calling z uzyciem freebayes i bcftools

find Projekt/bams -name "*.bam" > Projekt/bams/bams_list.txt

freebayes -f Projekt/ref/R64-1-1.fasta \
    -L Projekt/bams/bams_list.txt \
    > Projekt/freebayes/freebayes_variants.vcf

bcftools mpileup -Ou -f Projekt/ref/R64-1-1.fasta \
    --bam-list Projekt/bams/bams_list.txt \
    | bcftools call -mv -Ov -o Projekt/bcftools/bcftools_variants.vcf



