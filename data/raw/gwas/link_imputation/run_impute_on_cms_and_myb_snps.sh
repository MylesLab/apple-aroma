#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=80G
#SBATCH --account=def-smyles
#SBATCH --job-name=cms-myb-link-impute
#SBATCH --output=cms-myb-link-impute.out
#SBATCH --time=14-1:00:00

#SBATCH --mail-user=tayab.soomro@dal.ca
#SBATCH --mail-type=ALL

export JAVA_TOOL_OPTIONS="-Xmx60g"
export _JAVA_OPTIONS="-Xmx60g"

java -Xmx60g -jar /project/def-smyles/myles_lab/bin/LinkImpute/latest/LinkImpute.jar -q /project/def-smyles/myles_lab/2021-GCMS-Project/analysis/adding_snps_to_gbs/Markers/cms_myb_markers.ped ~/scratch/adding_snps_to_gbs/link_imputation/cms_myb_markers.imputed.ped
