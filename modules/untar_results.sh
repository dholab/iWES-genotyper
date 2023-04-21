#!/usr/bin/env bash

echo $1

echo $@

for arg in "$@"; do
  shift
  case "$arg" in
    "--out_dir") set -- "$@" "-o" ;;
    "--in_dir") set -- "$@" "-i" ;;
    "--submit_name") set -- "$@" "-s" ;;
    "--script_dir") set -- "$@" "-d" ;;
    *)        set -- "$@" "$arg"
  esac
done

while getopts o:i:s:d: opt ; do
   case $opt in
	    o) out_dir=$OPTARG ;;
	    i) in_dir=$OPTARG ;;
      s) submit_name=$OPTARG ;;
      d) script_dir=$OPTARG ;;
      *) usage; exit 1;;
   esac
done




find ${in_dir} -maxdepth 2 -type f -name "*.tar" > ${in_dir}/result_tar.txt
cd ${out_dir}
for tar_file in `cat ${in_dir}/result_tar.txt`; do
  echo ${tar_file}
  tar -xvf ${tar_file}
done

# concat result files
python3 ${script_dir}/concat_result_files.py --in_dir=${out_dir} --submit_name=${submit_name}