#!/usr/bin/env bash

usage () { echo "Usage : $0 "; }
echo $1

echo $@
create_reference_files=T
threads=1
ram=8000
script_dir=NONE
read_length=151
for arg in "$@"; do
  shift
  case "$arg" in
    "--ref_dir") set -- "$@" "-c" ;;
    "--cp_dir") set -- "$@" "-p" ;;
    "--fastq_dir") set -- "$@" "-f" ;;
    "--threads") set -- "$@" "-t" ;;
    "--ram") set -- "$@" "-r" ;;
    "--out_dir") set -- "$@" "-o" ;;
    "--project_name") set -- "$@" "-n" ;;
    "--animal_lookup_path") set -- "$@" "-a" ;;
    "--create_reference_files") set -- "$@" "-i" ;;
    "--db_path") set -- "$@" "-d" ;;
    "--db_ext") set -- "$@" "-e" ;;
    "--haplo_fasta") set -- "$@" "-m" ;;
    "--haplotype_json_path") set -- "$@" "-h" ;;
    "--script_dir") set -- "$@" "-z" ;;
    "--species") set -- "$@" "-s" ;;
    "--use_mask") set -- "$@" "-u" ;;
    "--output_depth_all") set -- "$@" "-b" ;;
    "--confidence_coeff_json") set -- "$@" "-g" ;;
    "--read_length") set -- "$@" "-l" ;;
    "--minimap2_path") set -- "$@" "-y" ;;
    *)        set -- "$@" "$arg"
  esac
done


while getopts c:p:f:t:r:o:n:a:i:d:e:m:h:z:s:u:b:g:l:y: opt ; do
   case $opt in
      c) ref_dir=$OPTARG ;;
	    p) cp_dir=$OPTARG ;;
	    f) fastq_dir=$OPTARG ;;
	    t) threads=$OPTARG ;;
      r) ram=$OPTARG ;;
	    o) out_dir=$OPTARG ;;
	    n) project_name=$OPTARG ;;
      a) animal_lookup_path=$OPTARG ;;
      i) create_reference_files=$OPTARG ;;
      d) db_path=$OPTARG ;;
      e) db_ext=$OPTARG ;;
      m) haplo_fasta=$OPTARG ;;
      h) haplotype_json_path=$OPTARG ;;
      z) script_dir=$OPTARG ;;
      s) species=$OPTARG ;;
      u) use_mask=$OPTARG ;;
      b) output_depth_all=$OPTARG ;;
      g) confidence_coeff_json=$OPTARG ;;
      l) read_length=$OPTARG ;;
      y) minimap2_path=$OPTARG ;;
      *) usage; exit 1;;
   esac
done

if [ "$script_dir" == "NONE" ]; then
  script_dir=$( dirname -- "$0"; );
fi
echo "Script dir: ${script_dir}"
# Launch the ipd ref matrix creator
create_reference_files=`echo $create_reference_files | cut -c1-1`
case "${create_reference_files}" in ([Tt])
  echo "python3 ${script_dir}/modules/create_ref_files.py \
--ref_dir=${ref_dir} \
--db_path ${db_path} \
--db_ext ${db_ext} \
--haplo_fasta=${haplo_fasta} \
--haplotype_json_path=${haplotype_json_path} \
--cp_path=${cp_dir} \
--threads=${threads} \
--ram=${ram} \
--species=${species} \
--minimap2_path=${minimap2_path}"
  python3 ${script_dir}/modules/create_ref_files.py \
  --ref_dir="${ref_dir}" \
  --db_path ${db_path} \
  --db_ext ${db_ext} \
  --haplo_fasta="${haplo_fasta}" \
  --haplotype_json_path="${haplotype_json_path}" \
  --cp_path="${cp_dir}" \
  --threads=${threads} \
  --ram=${ram} \
  --species=${species} \
  --minimap2_path=${minimap2_path};

  echo "python3 ${script_dir}/modules/create_hash.py --ref_dir=${ref_dir}";
  python3 ${script_dir}/modules/create_hash.py --ref_dir="${ref_dir}" --read_length=${read_length};
esac

## Align files
echo "python3 ${script_dir}/modules/semiperfect_align.py \
--cp_dir=${cp_dir} \
--fastq_dir=${fastq_dir} \
--bam_dir=${out_dir}/bam \
--ref_dir=${ref_dir} \
--threads=${threads} \
--ram=${ram}"
python3 ${script_dir}/modules/semiperfect_align.py \
--cp_dir="${cp_dir}" \
--fastq_dir="${fastq_dir}" \
--bam_dir="${out_dir}/bam" \
--ref_dir="${ref_dir}" \
--threads=${threads} \
--ram=${ram}

# Expand the alignment and create depth of coverage plots
echo "python3 ${script_dir}/modules/filter_alignments.py \
--project_name=${project_name} \
--out_dir=${out_dir} \
--bait_fasta=${ref_dir}/bait.fasta \
--ipd_ref_hash=${ref_dir}/ipd_hash.csv.gz \
--bam_dir=${out_dir}/bam \
--mask_path=${ref_dir}/mask_range.csv \
--db_ext ${db_ext} \
--output_depth_all ${output_depth_all} \
--use_mask ${use_mask}"


python3 ${script_dir}/modules/filter_alignments.py \
--project_name=${project_name} \
--out_dir=${out_dir} \
--bait_fasta=${ref_dir}/bait.fasta \
--ipd_ref_hash=${ref_dir}/ipd_hash.csv.gz \
--bam_dir=${out_dir}/bam \
--mask_path=${ref_dir}/mask_range.csv \
--db_ext ${db_ext} \
--output_depth_all ${output_depth_all} \
--use_mask ${use_mask}

# create a pivot table
echo "python3 ${script_dir}/modules/create_pivot.py \
      --project_name=${project_name} \
      --out_dir=${out_dir} \
      --ref_dir=${ref_dir} \
      --confidence_coeff_json=${confidence_coeff_json} \
      --animal_lookup_path=${animal_lookup_path} \
      --db_ext ${db_ext}"

python3 ${script_dir}/modules/create_pivot.py \
--project_name="${project_name}" \
--out_dir="${out_dir}" \
--ref_dir="${ref_dir}" \
--confidence_coeff_json=${confidence_coeff_json} \
--animal_lookup_path="${animal_lookup_path}" \
--db_ext ${db_ext}