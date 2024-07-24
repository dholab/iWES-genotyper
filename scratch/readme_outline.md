# IWES Genotyping 
## Purpose
- This pipeline is for the alignment and calling of MHC Class I and Class II alleles for non-human primates using high-accuracy short-length paired end reads
  - Currently optimized for Mamu, Mafa
  - Uses Curated reference databases from IPD
  - Currently optimized for Illumina NovaSeq reads with Twist or SeqCap Capture probes with reads paired end 151 BP
    - Can be used with longer reads (i.e. 251 bp) but is untested
  - Can be used for other regions but would need optimizations.
## Challenges
- The MHC in non-human primates is highly recombinant, and can have very similar alleles (as little as 1 snp difference)
- Each animal can have dozens of similar alleles
  - Many reads can map to multiple animals
- Many alleles are chimeras of 1 or more different real alleles. (A deletion or insertion in a repeat bp sequence would be chimeric)
  - Chimeric alignment alleles can have full depth of coverage despite not being "really" in the sample 
    - This is extremely common, nearly every animal has at least one chimeric allele
    - This full coverage comes from reads generated from other alleles in the sample, that when aligned to a Chimeric alignment allele can fully cover the allele
- The depth of coverage for Twist Probes is consistent from sample to sample allele type to allele type.
- Close to the edges have reduced depth of coverage due to semiperfect alignment criteria
- Valid alleles may exist in a sample that are not part of the reference database
## Approaches
- Retain all reads that semiperfectly align (even if they semiperfectly align to multiple alleles) 
  - Expand all semiperfect alignments to every allele, as long as the semiperfect criteria is met
    - This means some reads can have a reduced score (150 of 151 bp mapped on one but 151 of the other because the edge of the reference sequence is different)
    - But as long as they meet the semiperfect criteria, not just score, they should be retained.
- Ensure that the Reverse Paired end also maps unless:
  - The Read is close to the edge (within 450 bp) AND no other allele has both the forward and reverse paired end mapped.
- Allow for zero depth of coverage close to the edges (within ~40 basepairs to the edge)
  - Document for each read how much of the edges were excluded.
- Determine the maximum gap between consecutive alignment positions for mapped alleles (based on start or end positions).
  - Ideally, the start position of each consecutive already be 1, meaning every start and end position would be covered
  - However that is statistically, astronomically low probability of occurring with our protocol (<1E-30) 
  - If this gap is large, it may signify a chimeric alignment.
  - Therefore, filter if the gap is more that 80 bp between positions of consecutive alignments
- Document Passing/Called alleles:
  - Alleles must have a depth of coverage of at least 3 reads in every position of except for the excluded regions (i.e. close to the edge)
- Track how many "passing" alleles each read (or read pair) semiperfectly aligns to.
- Track if the allele has at least one read that uniquely aligns to the allele, and no other "passing" allele
- Normalize the read count to number of "passing" alleles the reads align to. (Contributes 1/N mappings: 1 for one allele, 1/2 for two alleles etc.)
- Track the Max gap for each allele
- Track the positions from the each edge that is excluded from the depth of coverage requirement
- Track the normalize median depth of coverage for each allele across each position with at least 1 read aligned (exclude zero).
- Present results as a pivot table for additional analaysis and manual curation
- Generate BAM files and depth of coeverage files for the passing alleles to help with manual curation.
- Highlight samples with medians that are lower than expected (0.75 of the median)
- Highlight samples with depth of coverage edge exclusions
- Highlight samples with large max gaps
- Highlight samples without a uniquely mapped read
- In some cases, the exon region of a samples' allele may be in a reference database but not the entire genomic region.
  - For this purpose we allow multiple database reference fasta's to be used (with calling priority, most specific to least specific i.e. genomic to exon2 if the exon region is captured in the genomic region of a read that is also called)
- Create a script to quickly expand ALL alignments 
  - We used a very fast hash table that can expand to all semiperfect alignments to all other semiperfect in seconds
- Create a script to filter the bam file to our additional requirements
  - We were unable to find a program that successfully found and retained all alignments.
  - We were unable to find a program (out of the box) to filter to our needs with our atypical filtering requirements.
  - bbmap_skimmer, found multiple alignments, 
    - But would only report the highest score, (not all valid semiperfect alignments)
    - Example: 
      - Allele 1: AGTAGTAGAGATCCATC 
      - Allele 2 (Del-A at Position1, C16G substiution): GTAGTAGAGATCCATG
      - READ: AGTAGTAGAGATCCAT would report mapping to Allele 1 but not Allele 2. Allele 1 would have a higher score for more aligned BP.

# Running the pipeline
## Create a docker image:
### Pre requisites
These packages were downloaded manually for the docker image.
- download BBMAP_39.01.tar.gz and place in directory of this Dockerfile
  - https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download
- download minimap2-2.24_x64-linux.tar.bz2 and place in directory of this Dockerfile
  - https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
## Build the image
- First ensure docker is installed on your machine. You may need to use "sudo" commands if you do not have root privileges
```
cd ~/
git clone https://github.com/dholab/iWES-genotyper
cd ~/iWES-genotyper/docker
docker build -t iwes_genotyper:v1_0 .
# Assumes data is in /Volumes/T8/iwes and iWES-genotyper repository code is in cd /home/username/iWES-genotyper
docker run -it \
-v /Volumes/T8/iwes:/Volumes/T8/iwes \
-v /home/username/iWES-genotyper:/home/username/iWES-genotyper \
iwes_genotyper:v1_0

# now run your commands bellow
```

## Running the pipeline locally
- change the flags/variable declarations to point to your paths and set up on your computer
```shell
run_dir=/Volumes/T8/iwes/run_1
script_dir=~/iwes_genotyping_v2
input_ref_dir=${script_dir}/ref/Mamu
ref_dir=/Volumes/T8/iwes/Mamu
minimap2_path=minimap2
cp_dir=~/anaconda3/bin/bbmap/current
cp_dir=/bbmap/current
${script_dir}/iwes_pipeline.sh \
  --ref_dir ${ref_dir} \
  --cp_dir ${cp_dir} \
  --fastq_dir ${run_dir}/fastq \
  --threads 4 \
  --ram 8000 \
  --out_dir ${run_dir}/out \
  --project_name 103126 \
  --animal_lookup_path ${run_dir}/run_1_animal_lookup.csv \
  --create_reference_files True \
  --db_path "${input_ref_dir}/Mamu_MHC-genomic.fasta ${input_ref_dir}/Mamu_MHC-exon2.fasta ${input_ref_dir}/Mamu_MHC-miseq.fasta" \
  --db_ext "ipd exon miseq" \
  --haplo_fasta ${input_ref_dir}/26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta \
  --haplotype_json_path ${input_ref_dir}/haplotype_lookup.json \
  --script_dir ${script_dir} \
  --species Mamu \
  --use_mask False \
  --output_depth_all False \
  --confidence_coeff_json ${input_ref_dir}/confidence_coeff.json \
  --minimap2_path ${minimap2_path}
```
## Flags: 
- minimap2_path (str): optional default: "minimap2" | path where your call minimap2.
- project_name (str): it is used as a prefix for files and folders.
- fastq_dir (str): mounted directory where your reverse paired end reads fastq.gz files reside.
- bam_dir (str): optional (default = ${out_dir}/bam ) | directory where you first alignment files will are outputted (semi perfect mode to your bait.fasta)
- out_dir (str): directory where your aggregated and intermediate files are outputted.
- animal_lookup_path (optional): path to the animal lookup csv file with headers and first column animal_id and second column your gs_id (file number) prefix of the paired fastq files
- ref_dir (str): directory of your aggregated reference files.
- db_path (str): Space delimited list of reference database fasta files. In order of prioirty (most specific i.e. genomic to least specific miseq).
- db_ext (str): Space delimited list reference database extensions i.e. "gen exon miseq" or "gen" or "gen exon".
- haplo_fasta (str): optional | filepath to your the reference database fasta your haplotype classifying json uses (fasta).
- haplotype_json_path: optional | filepath to a curated json dictionary to convert from called alleles to classifying haplotypes.
- cp_dir (str): file path to where bbmap align2 file is stored depends on the install method of bbmap.
- threads (int):default: 1 | number of threads (cpu's) to allocate to the alignment step
- ram (int): default: 8000 | amount of ram (in MB) to allocate to the alignment step
- script_dir (str): path of your iwes_genotyping_v2 directory
- species (str): must match prefix in your reference files (i.e. Mamu, Mafa)
- read_length (int): optional default: 151| the minimum length of your reads you want to use (often length of read)
- create_reference_files (bool: True, False, T, F, t,f): default: True | if you ran it with T prior and the files \
and is sotred in your config_dir, then you can set it to false. \
you need to run it for each set of reference files (i.e. different species, updated reference fasta databases) \
it must be reran even if only character line of one reference file is changed
- output_depth_all (bool: True, False, T, F, t,f): default: False| Create depth of coverage for all alleles with at least 1 aligned read, not just passing alleles
- use_mask (bool): default: True| file in path ${ref_dir}/mask_range.csv. By default the mask_range.csv generated is blank and will change the filter
# Running pipeline with serpant in Center for High-Throughput Computing - CHTC
- You must have access to Wisc.edu CHTC.  
  - If you do not have access this will not work and use the result-equivalent local pipeline
  - The only advantage of using CHTC is faster turn around of results by using a distributive computing
  - The final results should be identical (intermediate results generate additional files required for distributed computing)
- Get the latest version of serpent:
- https://github.com/DABAKER165/chtc_serpent_v2
## Create default and chtc configuration files for your runs
- Change the settings for the variables at the top of the script
```shell
###############################################
# These configure defaults that seldom change #
serpent_script_dir=~/chtc_serpent_v2
local_serpent_home_parent_dir=~/
input_run_dir=~/input_serpent
local_serpent_in_parent_dir=/Volumes/T8
local_serpent_out_parent_dir=${local_serpent_in_parent_dir}
submit_username=username
infiles_server_username=${submit_username}
in_files_server=yourtransferserver.domain.edu
out_files_server=${in_files_server}
in_file_home_dir=/staging/your/path
submit_server=yoursubmitserver.domain.edu
submit_home_dir=/home/${submit_username}

# More CHTC Node (submit file) Settings
docker_server=your.domain.edu
docker_user=username
docker_image=iwes
docker_tag=27329
docker_image_string=${docker_server}/${docker_user}/${docker_image}:${docker_tag}
chtc_ram=16GB
chtc_cpus=1
chtc_disk_space=30GB
# machine_requirements="Target.HasCHTCStaging == true"
# add any other flags (i.e. if you have your own priority server)
# "priority_flag": "accounting_group = yourgroup",
machine_requirements=your_requirement_string
# End configure defaults that seldom change #
#############################################
#########################################
# Start Submit file/ Run based settings #
minimap2_path=minimap2
#cp dir for server and/or docker image for server:
cp_path_server=/opt/conda/opt/bbmap-38.90-3/current
#cp dir for local and/or docker image for local:
cp_path=~/anaconda3/bin/bbmap/current

project_name=run1 #same as project_name
skip_create_ref_files=False
# assumes you saved it to your home directory i.e. /home/username/ on debian/ubuntu or similar; /Users/username/ on MacOSx
local_pipeline_code_dir=~/iWES-genotyper
default_config_dir=${local_pipeline_code_dir}/config_mamu_class1_2
species=Mamu
input_ref_dir=${local_pipeline_code_dir}/ref/Mamu
ref_dir=/Volumes/T8/iwes/Mamu
db_path="${input_ref_dir}/Mamu_MHC-genomic.fasta ${input_ref_dir}/Mamu_MHC-exon2.fasta ${input_ref_dir}/Mamu_MHC-miseq.fasta" 
db_ext="gen exon miseq"
haplo_fasta=${input_ref_dir}/Mamu_MHC-miseq.fasta
haplotype_json_path=${input_ref_dir}/haplotype_lookup.json

threads=4 #(4 vcpus, do not over provision)
ram=8000 # (8000 MB)
run_dir=/Volumes/T8/iwes/run_1
input_fastq_dir=${run_dir}/fastq
animal_lookup_path=${run_dir}/run_1_animal_lookup.csv
read_length=151
output_depth_all=False
use_mask=False
output_depth_all=False

# End Submit file/ Run based settings #
#######################################
############################################################################################################
### Do not change below this line (unless attempting develop settings not described in this repository #####

mkdir -p ${default_config_dir}
mkdir -p ${input_run_dir}/config/

echo "{
  \"all\": {
    \"local\": {
      \"submit_paths\": {
        \"home_dir\": \"${local_serpent_home_parent_dir}\",
        \"pipeline_code_dir\": \"${local_pipeline_code_dir}\"
      },
      \"in_paths\": {
        \"home_dir\": \"${local_serpent_in_parent_dir}\"
      },
      \"out_paths\": {
        \"home_dir\": \"${local_serpent_out_parent_dir}\"
      }
    },
    \"chtc\": {
      \"submit_paths\": {
        \"un\": \"${submit_username}\",
        \"server\": \"${submit_server}\",
        \"home_dir\": \"${submit_home_dir}\"
      },
      \"in_paths\": {
        \"un\": \"${infiles_server_username}\",
        \"server\": \"${in_files_server}\",
        \"home_dir\": \"${in_file_home_path}\"
      }
    }
  },
  \"create_ref_files\": {
    \"local\": {
      \"mark_as_completed\": \"${skip_create_ref_files}\",
      \"executable\": \"create_ref_files.py\",
      \"arguments\": {
        \"--ref_dir\": \"${ref_dir}\",
        \"-d\": \"${db_path}\",
        \"--db_ext\": \"${db_ext}\",
        \"--haplo_fasta\": \"${haplo_fasta}\",
        \"--haplotype_json_path\": \"${haplotype_json_path}\",
        \"--cp_path\": \"${cp_path}\",
        \"--threads\": \"${threads}\",
        \"--species\": \"${species}\",
        \"--ram\": \"${ram}\",
        \"--minimap2_path\": \"${minimap2_path}\"
      }
    }
  },
  \"create_hash\": {
    \"local\": {
      \"mark_as_completed\": \"<create_ref_files:local:mark_as_completed>\",
      \"start_trigger\": \"<create_ref_files:completed>\",
      \"executable\": \"create_hash.py\",
      \"arguments\": {
        \"--ref_dir\": \"<create_ref_files:local:arguments:--ref_dir>\",
        \"--read_length\": \"${read_length}\"
      }
    }
  },
  \"iwes\": {
    \"local\": {
      \"start_trigger\": \"<create_hash:completed>\",
      \"mark_as_completed\": \"False\"
    },
    \"chtc\": {
      \"submit_job\": \"True\",
      \"get_output\": \"True\",
      \"executable\": \"align_filter_chtc.sh\",
      \"transfer_to_server\": \"${input_fastq_dir}\",
      \"tar_files\": \"False\",
      \"sample_extension\": [
        \"_R1_001.fastq.gz\",
        \"_R2_001.fastq.gz\"
      ],
      \"static_files\": {
        \"ipd_ref_hash\": [
          \"<create_ref_files:local:arguments:--ref_dir>\",
          \"ipd_hash.csv.gz\"
        ],
        \"ref_fasta\": [
          \"<create_ref_files:local:arguments:--ref_dir>\",
          \"bait.fasta\"
        ],
        \"main_exon\": [
          \"<iwes:local:submit_paths:module_code_dir>\",
          \"filter_alignments.py\"
        ],
        \"mask_path\": [
          \"<create_ref_files:local:arguments:--ref_dir>\",
          \"mask_range.csv\"
        ]
      },
      \"arguments\": {
        \"--cp_dir\": \"${cp_path_server}\",
        \"--threads\": \"<iwes:chtc:cpus>\",
        \"--ram\": \"<iwes:chtc:ram>\",
        \"--project_name\": \"<submit_name>\",
        \"--db_ext\": \"<create_ref_files:local:arguments:--db_ext>\",
        \"--output_depth_all\": \"${output_depth_all}\"
      }
    }
  },
  \"untar_results\": {
    \"local\": {
      \"mark_as_completed\": \"False\",
      \"executable\": \"untar_results.sh\",
      \"start_trigger\": \"<iwes:completed>\",
      \"arguments\": {
        \"--out_dir\": \"<create_pivot:local:module_out_dir>\",
        \"--in_dir\": \"<iwes:local:module_out_dir>\",
        \"--submit_name\": \"<submit_name>\",
        \"--script_dir\": \"<untar_results:local:submit_paths:module_code_dir>\"
      }
    }
  },
  \"create_pivot\": {
    \"local\": {
      \"mark_as_completed\": \"False\",
      \"executable\": \"create_pivot.py\",
      \"start_trigger\": \"<untar_results:completed>\",
      \"arguments\": {
        \"--out_dir\": \"<create_pivot:local:module_out_dir>\",
        \"--ref_dir\": \"<create_ref_files:local:arguments:--ref_dir>\",
        \"--animal_lookup_path\": \"${animal_lookup_path}\",
        \"--project_name\": \"<submit_name>\",
        \"--db_ext\": \"<create_ref_files:local:arguments:--db_ext>\"
      }
    }
  }
}"> ${default_config_dir}/default.json

echo "{
  \"iwes\": {
    \"chtc\": {
      \"docker_image\": \"${docker_image_string}\",
      \"ram\": \"${chtc_ram}\",
      \"cpus\": \"${chtc_cpus}\",
      \"machine_requirements\": \"${machine_requirements}\",
      \"disk_space\": \"${chtc_disk_space}\"
    }
  }
}" > ${default_config_dir}/chtc.json

```
## create run based json
```shell
#######################################
# Start Submit file/ Run based settings #
skip_create_ref_files=False
local_pipeline_code_dir=~/iwes_genotyping_v2
default_config_dir=${local_pipeline_code_dir}/config_mamu_class1_2
species=Mamu
input_ref_dir=${local_pipeline_code_dir}/ref/Mamu
ref_dir=/Volumes/T8/iwes/Mamu
db_path="${input_ref_dir}/gen_exon_rem_2023_03_06.fasta ${input_ref_dir}/Mamu_MHC-II_Ex2_447_Seq_renamed_6Mar23.fasta ${input_ref_dir}/26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta" 
db_ext="gen exon miseq"
haplo_fasta=${input_ref_dir}/26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta
haplotype_json_path=${input_ref_dir}/haplotype_lookup.json
cp_path=~/anaconda3/bin/bbmap/current
threads=4 #(4 vcpus, do not over provision)
ram=8000 # (8000 MB)
run_dir=/Volumes/T8/iwes/run_1
input_fastq_dir=${run_dir}/fastq
animal_lookup_path=${run_dir}/run_1_animal_lookup.csv
cp_path_server=/opt/conda/opt/bbmap-38.90-3/current
read_length=151
output_depth_all=False
use_mask=False
# End Submit file/ Run based settings #
#######################################
############################################################################################################
### Do not change below this line (unless attempting develop settings not described in this repository #####
echo "{
  \"all\": {
    \"local\": {
      \"default_config_dir\": \"${default_config_dir}\"
    }
  },
    \"create_ref_files\": {
    \"local\": {
      \"mark_as_completed\": \"${skip_create_ref_files}\",
      \"arguments\": {
        \"--ref_dir\": \"${ref_dir}\",
        \"-d\": \"${db_path}\",
        \"--db_ext\": \"${db_ext}\",
        \"--haplo_fasta\": \"${haplo_fasta}\",
        \"--haplotype_json_path\": \"${haplotype_json_path}\",
        \"--cp_path\": \"${cp_path}\",
        \"--threads\": \"${threads}\",
        \"--species\": \"${species}\",
        \"--ram\": \"${ram}\",
        \"--minimap2_path\": \"${minimap2_path}\"
      }
    }
  },
  \"iwes\": {
	    \"local\": {
      \"mark_as_completed\": \"False\"
    },
    \"chtc\": {
      \"transfer_to_server\": \"${input_fastq_dir}\",
      \"arguments\": {
        \"--cp_dir\": \"${cp_path_server}\",
        \"--output_depth_all\": \"${output_depth_all}\",
        \"--output_depth_all\": \"${use_mask}\"
      }
    }
  },
  \"create_pivot\": {
    \"local\": {
    \"mark_as_completed\": \"False\",
    \"arguments\": {
        \"--animal_lookup_path\": \"${animal_lookup_path}\"
        }
    }
  }
}">${input_run_dir}/config/${project_name}.json

```
## Run chtc serpent
```shell
# modify chtc for priority flag as needed
# serpent_script_dir=~/chtc_serpent_v2
# generate your python run script
echo "python3 ${serpent_script_dir}/main.py \
--config_dir ${input_run_dir}/config \
--submission_dir ${input_run_dir}/chtc_serpent_submissions"
# copy and paste run script
```

# Create a mask ranges to exclude ranges of positions for each allele from the depth of coverage criteria:
## Purpose:
- Some alleles do not get complete coverage with a given capture probe system.
  - This could be because the probes were not included, or simply did not capture well
- This section allows you to maks certain ranges in an effort to maximize the amount of the genome used and increase the specificity of the calls
  Added the ability to crate a mask_range positions to exclude from the depth of coverage for each allele
## Creating a file manually:
- You can custom make a range mask, which are regions that will be exluded from depth of coverage requirements
- It is a 3 column .csv (comma-sperated values) file called "mask_range.csv" in the ref_dir
- Headers are case-sensitive: ALLELE, START, END
  - ALLELE is the exact header name in the fasta file of the bait.fasta reference db file.
  - START is the INCLUSIVE start position of the range (i.e. 1 would exclude the first position 
  - END is the INCLUSIVE end position of the range (i.e. 1 would exclude the first position)
  - Example
    - This exludes the first position of the depth of coverage requirements N1N7_NKG2A-gen
      - N1N7_NKG2A-gen, 1, 1
    - This excludes positions 300-400 of the depth of coverage requirements N1N7_NKG2A-gen
      - N1N7_NKG2A-gen, 300, 400
      
## Dynamically creating the ranges from existing data:
- "-o": (str) required: space delimited string, each is the outputted folder of results
  - It should be the parent folder of "depth_all" generated by running the align_filter_chtc.sh or filter_alignments.py with "output_depth_all=True"
- "-d": (str) required: the db suffix you want to mask.
- "--result_dir": (str) required: where the intermediate result files are stored
- "--ref_dir": (str) required: where the mask_range.csv is stored as well as the nmask.fasta
  - This should be the same ref_dir used in align_filter_chtc.sh or filter_alignments.py
- "--prefix": (str) optional, default='': an optional prefix to add to the intermediate files
- "--threshold": (int) optional, default=27, the depth of coverage criteria to use based on the maximum position for each allele.
  - This should be higher than the depth of coverage threshold in align_filter_chtc.sh or filter_alignments.py (default=3).
  - This should be higher than the depth of coverage threshold in align_filter_chtc.sh or filter_alignments.py (default=3).
  - 27 was chosen after trying other defaults and getting more and less over calls for our KIR and NKG2 sets in our publications (pending).  
  - You will need to adjust this value to determine what is best for your set.
```shell
# example of using two outputted runs
python3 ~/iwes_genotyping_v2/modules/create_mask_ranges.py \
-o /Volumes/T8/serpent/iwes_genotyping_v2/run_1/out_files/create_pivot_table \
/Volumes/T8/serpent/iwes_genotyping_v2/run_2/out_files/create_pivot_table \
-d gen \
--result_dir /Volumes/T8/serpent/iwes_genotyping_v2/run_1/out_files/create_mask \
--ref_dir /Volumes/T8/iwes/Mafa_nkg2_mask \
--threshold 27

# example of using one outputted runs
python3 ~/iwes_genotyping_v2/modules/create_mask_ranges.py \
-o /Volumes/T8/serpent/iwes_genotyping_v2/run_1/out_files/create_pivot_table \
-d gen \
--result_dir /Volumes/T8/serpent/iwes_genotyping_v2/run_1/out_files/create_mask \
--ref_dir /Volumes/T8/iwes/Mafa_nkg2_mask \
--threshold 27
```
## How to apply mask ranges
- First run the main  work flow with your set. Back up files as needed.
- Next run the above code based on where your files are located
- Make sure the generate mask range file (i.e. mask_range.csv) is in your ref_dir
- Finally, rerun the main workflow with the mask range file in your ref_dir and it should pick up automatically.

# Acknowledgments:
## Software:
minimap2:
- Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705

samtools:
- Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

BBMap:
- Bushnell, Brian. BBMap: A Fast, Accurate, Splice-Aware Aligner. United States: N. p., 2014.

muscle:
- Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

### Python packages
Pandas:
```
@software{reback2020pandas,
    author       = {The pandas development team},
    title        = {pandas-dev/pandas: Pandas},
    month        = feb,
    year         = 2020,
    publisher    = {Zenodo},
    version      = {latest},
    doi          = {10.5281/zenodo.3509134},
    url          = {https://doi.org/10.5281/zenodo.3509134}
}
@InProceedings{ mckinney-proc-scipy-2010,
  author    = { {W}es {M}c{K}inney },
  title     = { {D}ata {S}tructures for {S}tatistical {C}omputing in {P}ython },
  booktitle = { {P}roceedings of the 9th {P}ython in {S}cience {C}onference },
  pages     = { 56 - 61 },
  year      = { 2010 },
  editor    = { {S}t\'efan van der {W}alt and {J}arrod {M}illman },
  doi       = { 10.25080/Majora-92bf1922-00a }
}
```

## Development team:
Baker, David A <sup>1*</sup>, Minor, Nick R <sup>1</sup> Wiseman, Roger <sup>1</sup>, Karl, Julie <sup>1</sup>, Prall, Trent <sup>1</sup>, O'Connor David H. <sup>1</sup>
- more names to be added.
- <sup>*</sup> denotes most significant contributions
- <sup>1</sup> University of Wisconsin-Madison, School of Public Health and Medicine
- <sup>2</sup> Baylor, School of Medicine
