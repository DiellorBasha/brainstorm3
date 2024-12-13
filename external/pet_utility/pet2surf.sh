#!/bin/bash

# PET to Surface Mapping Script
# =============================
# Usage:
# ------
# ./pet2surf.sh -s subject_label -p pet_dir [-f fwhm] [-n threads]
#
# Arguments:
# -s subject_label : The label of the subject in BIDS format (e.g., sub-01)
# -p pet_dir       : Path to the PET directory containing the subject's data
# -f fwhm          : (Optional) Full Width at Half Maximum (FWHM) for point spread function. Default is 2.8.
# -n threads       : (Optional) Number of threads for parallel processing. Default is 1.
# -h               : Display help

# Function to display usage information
show_help() {
  echo "Usage: ./pet2surf.sh -s subject_label -p pet_dir [-f fwhm] [-n threads]"
  echo "  -s subject_label : Subject label in BIDS format (e.g., sub-01)"
  echo "  -p pet_dir       : Path to the PET directory containing the subject's data"
  echo "  -f fwhm          : FWHM for point spread function (default: 2.8)"
  echo "  -n threads       : Number of threads for parallel processing (default: 1)"
  echo "  -h               : Display this help message"
}

# Default parameters
FWHM=6
THREADS=1

# Parse command-line arguments
while getopts "s:p:f:n:h" opt; do
  case ${opt} in
    s) SUBJECT_LABEL=$OPTARG ;;
    p) PET_DIR=$OPTARG ;;
    f) FWHM=$OPTARG ;;
    n) THREADS=$OPTARG ;;
    h) show_help; exit 0 ;;
    *) show_help; exit 1 ;;
  esac
done

# Check required arguments
if [ -z "$SUBJECT_LABEL" ] || [ -z "$PET_DIR" ]; then
  echo "Error: Missing required arguments."
  show_help
  exit 1
fi

# Check if SUBJECTS_DIR environment variable is set
if [ -z "$SUBJECTS_DIR" ]; then
  echo "Error: The SUBJECTS_DIR environment variable is not set."
  exit 1
fi

# Define directories
FS_PET_DIR="$PET_DIR"
T1="$SUBJECTS_DIR/$SUBJECT_LABEL/mri/T1.mgz"
PET_FILE="$FS_PET_DIR/PET.nii"
TEMPLATE="$FS_PET_DIR/PET.mgz"
SURF="$SUBJECTS_DIR/$SUBJECT_LABEL/surf"

# Ensure required input files exist
if [ ! -f "$PET_FILE" ]; then
  echo "Error: PET file not found: $PET_FILE"
  exit 1
fi

echo "========Conversion==============+"
    
# Convert PET file to FreeSurfer format
mri_convert "$PET_FILE" "$TEMPLATE" || exit 1

echo "========Coregistration==========="
# Register the template to generate a PET2MRI.reg.lta file
 if [ ! -f "$FS_PET_DIR/PET2MRI.reg.lta" ]; then
  echo "Running mri_coreg to generate PET2MRI.reg.lta..."
  mri_coreg --s "$SUBJECT_LABEL" \
            --mov "$TEMPLATE" \
            --reg "$FS_PET_DIR/PET2MRI.reg.lta" \
            --threads $THREADS
else
  echo "Skipping mri_coreg: PET2MRI.reg.lta already exists."
fi
 
 if [ ! -f "$FS_PET_DIR/PET_REG.mgz" ]; then
  echo "Running mri_vol2vol to generate PET_REG.mgz..."
  mri_vol2vol --mov "$TEMPLATE" \
                --targ "$T1" \
                --lta "$FS_PET_DIR/PET2MRI.reg.lta" \
                --o "$FS_PET_DIR/PET_REG.mgz"
else
  echo "Skipping mri_vol2vol: PET_REG.mgz already exists."
fi

# Apply smoothing
  mri_convert --fwhm $FWHM \
              "$FS_PET_DIR/PET_REG.mgz" \
              "$FS_PET_DIR/PET_REG_SMOOTH.mgz" \
              --nthreads $THREADS


if [ ! -f "$FS_PET_DIR/PET_REG_MASKED.mgz" ]; then
  echo "Running mri_mask to generate PET_REG_MASKED.mgz..."
 # Mask the PET volume
  mri_mask "$FS_PET_DIR/PET_REG_SMOOTH.mgz" \
          "$SUBJECTS_DIR/$SUBJECT_LABEL/mri/brainmask.mgz" \
          "$FS_PET_DIR/PET_REG_MASKED.mgz" || exit 1
else
  echo "Skipping mri_mask: PET_REG_MASKED.mgz already exists."
fi
 

echo "========Mapping to Surface==========="

  mri_vol2surf --src "$FS_PET_DIR/PET_REG_MASKED.mgz" \
              --srcreg "$FS_PET_DIR/PET2MRI.reg.lta" \
              --hemi lh \
              --projfrac 1 \
              --o  "$FS_PET_DIR/lh_pet.nii.gz" \
              --no-reshape || exit 1

  mri_vol2surf --src "$FS_PET_DIR/PET_REG_MASKED.mgz" \
              --srcreg "$FS_PET_DIR/PET2MRI.reg.lta" \
              --hemi rh \
              --projfrac 1 \
              --o  "$FS_PET_DIR/rh_pet.nii.gz" \
              --no-reshape || exit 1

echo "========Converting Surfaces to Curv==========="

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval "$FS_PET_DIR/lh_pet.nii.gz" \
  --trgsubject "$SUBJECT_LABEL" \
  --trgsurfval "$FS_PET_DIR/petproj" \
  --hemi lh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval "$FS_PET_DIR/rh_pet.nii.gz" \
  --trgsubject "$SUBJECT_LABEL" \
  --trgsurfval "$FS_PET_DIR/petproj" \
  --hemi rh \
  --trg_type curv
echo "========Processing completed for subject $SUBJECT_LABEL==========="



