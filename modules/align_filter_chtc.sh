#!/usr/bin/env bash

echo $1

threads=1
ram=16000
ref_fasta=bait.fasta
db_ext="gen exon miseq"
ipd_ref_hash=ipd_hash.csv.gz
output_depth_all=False
mask_path=mask_range.csv
use_mask=True
for arg in "$@"; do
  shift
  case "$arg" in
    "--cp_dir") set -- "$@" "-p" ;;
    "--threads") set -- "$@" "-t" ;;
    "--ram") set -- "$@" "-r" ;;
    "--project_name") set -- "$@" "-n" ;;
    "--ipd_ref_hash") set -- "$@" "-f" ;;
    "--ref_fasta") set -- "$@" "-a" ;;
    "--db_ext") set -- "$@" "-x" ;;
    "--mask_path") set -- "$@" "-m" ;;
    "--use_mask") set -- "$@" "-u" ;;
    "--output_depth_all") set -- "$@" "-o" ;;
    *)        set -- "$@" "$arg"
  esac
done

while getopts p:t:r:b:n:f:a:x:m::o:u: opt ; do
   case $opt in
	    p) cp_dir=$OPTARG ;;
	    t) threads=$OPTARG ;;
      r) ram=$OPTARG ;;
	    n) project_name=$OPTARG ;;
      f) ipd_ref_hash=$OPTARG ;;
      a) ref_fasta=$OPTARG ;;
      x) db_ext=$OPTARG ;;
      m) mask_path=$OPTARG ;;
      u) use_mask=$OPTARG ;;
      o) output_depth_all=$OPTARG ;;
      *) usage; exit 1;;
   esac
done

# Ram value expected (int in mb, or ends with MB, GB, or TB.)
# get last 2
ram_unit=$(echo -n $ram | tail -c 2)
# the 2nd to last it the identifier
ram_unit=$(echo -n $ram_unit | cut -c1-1)

case "${ram_unit}" in
  [Tt])
    ram=${ram%??}
    ram=$(( 1000000 * ram ));;
  [Gg])
    ram=${ram%??}
    ram=$(( 1000 * ram ));;
  [Mm])
    ram=${ram%??};;
  *)
    echo -n "RAM variable bad format: RAM value expected: int in mb, or ends with MB, GB, or TB. i.e. \"8GB\")"
esac

#ref_tar_file=$( basename -- "$ref_tar"; )
#tar -xvf ${ref_tar_file}
#rm -f ${ref_tar_file}
#ref_tar_dir=${ref_tar_file%????}

for in_1 in *_R1_001.fastq.gz; do
  first_char=$(echo "${in_1}" | cut -c1-1);
  if [[ "$first_char" == "." ]]; then
    continue
  else
    break
  fi
done

file_prefix=${in_1%????????????????}
in_2=${file_prefix}_R2_001.fastq.gz


bam_dir=bam
mkdir -p "${bam_dir}"

java -ea -Xmx${ram}m -Xms${ram}m -cp ${cp_dir} align2.BBMap build=1 \
    in=${in_1} in2=${in_2} \
    ref=./${ref_fasta} \
    outm=${bam_dir}/${file_prefix}.bam \
    semiperfectmode=t \
    threads=${threads} \
    nodisk=f

echo "finished bbmap semiperfect alignment of sample: ${file_prefix} in project: ${project_name}"
    # '.format(int(ram), cp_dir, in_1, in_2, bait_fasta, bam_dir, sample_i, int(threads)))
# Expand the alignment and create depth of coverage plots
echo "python3 filter_alignments.py \
--project_name=${project_name} \
--out_dir=\".\" \
--ipd_ref_hash=\"./${ipd_ref_hash}\" \
--bam_dir=${bam_dir} \
--bait_fasta=./${ref_fasta} \
--mask_path=./${mask_path} \
--db_ext ${db_ext} \
--output_depth_all ${output_depth_all} \
--use_mask ${use_mask}"

python3 filter_alignments.py \
--project_name=${file_prefix} \
--out_dir=. \
--bait_fasta=./${ref_fasta} \
--ipd_ref_hash=./${ipd_ref_hash} \
--bam_dir=${bam_dir} \
--mask_path=./${mask_path} \
--db_ext ${db_ext} \
--output_depth_all ${output_depth_all} \
--use_mask ${use_mask}

rm -rf ipd_ref_matrix
rm -f ${ipd_ref_hash}
rm -f ${ref_fasta}
rm -f main_exon.py
rm -f ${ref_fasta}.fai
rm -f ${in_1}
rm -f ${in_2}

tar -cvf ${file_prefix}.tar \
${file_prefix}.depth_all.csv \
${file_prefix}_norm_median_all.csv \
${file_prefix}_read_ct.csv \
${file_prefix}_gaps.csv \
depth-exon \
depth-gen \
depth-miseq \
exon-bam \
bam_all-gen \
bam_all-exon \
bam_all-miseq \
bam-gen \
bam-exon \
bam-miseq \
bam \
expanded_maps-exon \
expanded_maps-gen \
expanded_maps-miseq \
normalized_median-exon \
normalized_median-gen \
normalized_median-miseq
rm -f ${file_prefix}_norm_median-exon.csv
rm -f ${file_prefix}_norm_median-gen.csv
rm -f ${file_prefix}_norm_median-miseq.csv
rm -f ${file_prefix}_norm_median_all.csv
rm -f ${file_prefix}_read_ct.csv
rm -f ${file_prefix}_gaps.csv
