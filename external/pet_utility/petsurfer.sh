#!/bin/bash

# PETSurfer Pipeline Script
# =========================
# Usage:
# ------
# ./petsurfer.sh -s subject_label -p pet_dir [-f fwhm] [-n threads]
#
# Arguments:
# -s subject_label : The label of the subject in BIDS format (e.g., sub-01)
# -p pet_dir       : Path to the BIDS PET directory containing the subject data
# -r ref_region    : Reference region for SUVR calculation (e.g., Cerebellum-Cortex)
# -f fwhm          : (Optional) Full Width at Half Maximum (FWHM) for PSF in PVC. Default is 2.8.
# -n threads       : (Optional) Number of threads for parallel processing. Default is 1.
# -h      

show_help() {
  echo "Usage: ./petsurfer.sh -s subject_label -p pet_dir [-f fwhm] [-n threads] [-r ref_region] [-x]"
  echo "  -s subject_label : Subject label in BIDS format (e.g., sub-01)"
  echo "  -p pet_dir       : Path to PET data directory"
  echo "  -r ref_region    : Reference region (default: Cerebellum-Cortex)"
  echo "  -f fwhm          : FWHM for PVC (default: 2.8)"
  echo "  -n threads       : Number of threads (default: 1)"
  echo "  -x               : Disable Partial Volume Correction (PVC)"
  echo "  -h               : Display help"
}

USE_PVC=true # Default to using PVC
# Default FWHM and threads
FWHM=2.8 # Default FWHM for the point spread function (PSF) in PVC
THREADS=1 # Default number of threads for parallel processing
REF_REGION="Cerebellum-Cortex" # Default reference region for SUVR calculation

# Parse command-line arguments
while getopts "s:p:r:f:n:xh" opt; do
  case ${opt} in
    s) SUBJECT_LABEL=$OPTARG ;;
    p) PET_DIR=$OPTARG ;;
    r) REF_REGION=$OPTARG ;;
    f) FWHM=$OPTARG ;;
    n) THREADS=$OPTARG ;;
    x) USE_PVC=false ;; # Disable PVC if -x is provided
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
PVC="$FS_PET_DIR/pvc"
NOPVC="$FS_PET_DIR/nopvc"
# Path to FreeSurferColorLUT.txt
LUT_FILE="$FREESURFER_HOME/FreeSurferColorLUT.txt"

# Create necessary directories
mkdir -p "$PVC"
mkdir -p "$NOPVC"

# Read FreeSurferColorLUT.txt and find atlas values for a given region
get_atlas_values() {
  local region=$1
  local left=$(grep -w "Left-$region" "$LUT_FILE" | awk '{print $1}')
  local right=$(grep -w "Right-$region" "$LUT_FILE" | awk '{print $1}')
  
  if [ -z "$left" ] || [ -z "$right" ]; then
    echo "Error: Atlas values for region $region not found in $LUT_FILE"
    exit 1
  fi

  echo "$left $right"
}

# Get the atlas values for the selected reference region, or use defaults
if [ "$REF_REGION" == "Cerebellum-Cortex" ]; then
  REF_REGION_VALUES="8 47"
else
  REF_REGION_VALUES=$(get_atlas_values "$REF_REGION")
fi


echo "========Conversion==============+"
    
      mri_convert "$PET_FILE" \
                "$TEMPLATE"

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
 # mri_convert --fwhm $FWHM \
 #             "$FS_PET_DIR/PET_REG.mgz" \
 #             "$FS_PET_DIR/PET_REG_SMOOTH.mgz" \
 #             --nthreads $THREADS

if [ ! -f "$FS_PET_DIR/PET_REG.mgz" ]; then
  echo "Running mri_mask to generate PET_REG_MASKED.mgz..."
 # Mask the PET volume
  mri_mask "$FS_PET_DIR/PET_REG_SMOOTH.mgz" \
          "$SUBJECTS_DIR/$SUBJECT_LABEL/mri/brainmask.mgz" \
          "$FS_PET_DIR/PET_REG_MASKED.mgz" || exit 1
else
  echo "Skipping mri_mask: PET_REG_MASKED.mgz already exists."
fi

echo "========GTM Segmentation==========="
# Create new segmentation using FreeSurfer

# Define the path to the gtmseg output file
  gtmseg_output="$SUBJECTS_DIR/$SUBJECT_LABEL/mri/gtmseg.mgz"

# Check if gtmseg.mgz already exists
  if [ -f "$gtmseg_output" ]; then
      echo "gtmseg.mgz already exists for subject $SUBJECT_LABEL. Skipping gtmseg processing."
  else
      # Run the gtmseg command if the file does not exist
      echo "Running gtmseg for subject $SUBJECT_LABEL..."
      gtmseg --s "$SUBJECT_LABEL" --xcerseg || exit 1
  fi
  
if [ "$USE_PVC" = true ]; then
echo "========Partial Volume Correction (PVC)==========="
#if [ ! -f "$PVC/aux/bbpet2anat.lta" ];
  # Perform Partial Volume Correction (PVC)
  mri_gtmpvc --i "$TEMPLATE" \
            --reg "$FS_PET_DIR/PET2MRI.reg.lta" \
            --psf $FWHM \
            --seg gtmseg.mgz \
            --default-seg-merge \
            --auto-mask 1 .05 \
            --mgx .25 \
            --rescale $REF_REGION_VALUES \
            --save-input \
            --o "$PVC" \
            --threads $THREADS || exit 1
#else
 # echo "Found reg file."
#fi
echo "========Mapping PVC to Surface==========="

  mri_vol2surf --src "$PVC/mgx.ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi lh \
              --projfrac 1 \
              --o  "$SURF/lh_pet_pvc.nii.gz" \
              --no-reshape || exit 1

  mri_vol2surf --src "$PVC/mgx.ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi rh \
              --projfrac 1 \
              --o  "$SURF/rh_pet_pvc.nii.gz" \
              --no-reshape || exit 1

echo "========Converting Surfaces to Curv==========="

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/lh_pet_pvc.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $FS_PET_DIR/pet_pvc \
  --hemi lh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/rh_pet_pvc.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $FS_PET_DIR/pet_pvc \
  --hemi rh \
  --trg_type curv
echo "========Processing completed for subject $SUBJECT_LABEL==========="
echo "Surface projections: $FS_PET_DIR/pet_pvc"
else

echo "Skipping Partial Volume Correction (PVC) as per user input."
echo "==Calculating SUVR without Partial Volume Correction (noPVC)=="
if [ ! -f "$PVC/aux/bbpet2anat.lta" ]; then
  # Do mri_gtmpvc to obtain registration file 
  mri_gtmpvc --i "$TEMPLATE" \
            --reg "$FS_PET_DIR/PET2MRI.reg.lta" \
            --psf $FWHM \
            --seg gtmseg.mgz \
            --default-seg-merge \
            --auto-mask 1 .05 \
            --mgx .25 \
            --rescale $REF_REGION_VALUES \
            --save-input \
            --o "$PVC" \
            --threads $THREADS || exit 1
else
  echo "Found reg file."
fi
  # Get cortical mask
  mri_binarize --i $SUBJECTS_DIR/$SUBJECT_LABEL/mri/aseg.mgz \
             --match 3 42 \
             --o $PVC/ctxmask.mgz || exit 1
  
  # Mask rescaled PET image
  mri_mask $PVC/input.rescaled.nii.gz \
             $PVC/ctxmask.mgz \
             $NOPVC/ctxgm.nii.gz || exit 1


echo "========Mapping NOPVC to Surface==========="

  mri_vol2surf --src "$NOPVC/ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi lh \
              --projfrac 1 \
              --o  "$SURF/lh_pet_nopvc.nii.gz" \
              --no-reshape || exit 1

  mri_vol2surf --src "$NOPVC/ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi rh \
              --projfrac 1 \
              --o  "$SURF/rh_pet_nopvc.nii.gz" \
              --no-reshape

echo "========Converting Surfaces to Curv==========="

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/lh_pet_nopvc.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $FS_PET_DIR/pet_nopvc \
  --hemi lh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/rh_pet_nopvc.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $FS_PET_DIR/pet_nopvc \
  --hemi rh \
  --trg_type curv
echo "========Processing completed for subject $SUBJECT_LABEL==========="
echo "Surface projections: $FS_PET_DIR/pet_nopvc "

fi

