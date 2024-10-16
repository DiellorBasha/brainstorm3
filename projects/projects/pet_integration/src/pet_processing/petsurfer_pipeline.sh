#!/bin/bash

# PETSurfer Pipeline Script
# =========================
# 
# This script automates the processing of PET (Positron Emission Tomography) data using FreeSurfer's PETSurfer tools. 
# It processes PET data for a specific subject, aligning each PET frame to the subject's T1-weighted MRI, performing 
# partial volume correction (PVC), and mapping the PET data onto the cortical surface.
#
# Usage:
# ------
# ./petsurfer.sh -s subject_label -p pet_dir -t tracer [-f fwhm] [-n threads]
#
# Arguments:
# -s subject_label : The label of the subject in BIDS format (e.g., sub-01)
# -p pet_dir       : Path to the BIDS PET directory containing the subject data
# -t tracer        : Tracer name to match in the PET file (e.g., 18Fflortaucipir)
# -r ref_region    : Reference region for SUVR calculation (e.g., Cerebellum-Cortex)
# -f fwhm          : (Optional) Full Width at Half Maximum (FWHM) for PSF in PVC. Default is 2.8.
# -n threads       : (Optional) Number of threads for parallel processing. Default is 1.
# -h      
#
# Requirements and Assumptions:
# -----------------------------
# 1. The `SUBJECTS_DIR` environment variable must be set, pointing to the FreeSurfer subjects directory. 
#    SUBJECTS_DIR is typically set in the FreeSurfer environment setup script (e.g., freesurfer/SetUpFreeSurfer.sh)
#    and it would have been already set up for the reconstruction with recon-all. This is 
#    where the processed PET data and results will be stored and where the results of recon-all would have been stored. 
#    If in doubt, define your SUBJECTS_DIR in your terminal using
#         export SUBJECTS_DIR="/path/to/your/subjects_dir" 
#    and please check the FreeSurfer documentation.
# 2. The PET directory should follow the BIDS structure with files named according to the subject, session, 
#    and tracer (e.g., sub-01_ses-01_trc-18Fflortaucipir_pet.nii.gz).
# 3. The script assumes that FreeSurfer's recon-all has already been run for the subject, producing the necessary 
#    anatomical data in the `SUBJECTS_DIR`.
# 4. The PET data is registered to the subject's T1 MRI, and PVC is performed. The results are then mapped onto the 
#    cortical surface for further analysis.
# 5. Subject folders in the `SUBJECTS_DIR` should be named according to the BIDS format (e.g., sub-01) and should
#    contain the necessary anatomical data from FreeSurfer's recon-all. The subject directory must have this structure
#
#       |sub-01
#           |---mri (contains T1.mgz, brainmask.mgz, etc.)
#           |---surf (contains the surface files)
#           |---pet (will contain the processed PET data)
# 
# 6. The script uses FreeSurfer's tools (mri_convert, mri_coreg, mri_vol2vol, mri_gtmpvc, mri_vol2surf, mri_binarize,
#    mri_mask, gtmseg, mri_surf2surf) to process the PET data and generate the results.
# 7. The script assumes that the necessary FreeSurfer environment variables are set (e.g., $FREESURFER_HOME, $SUBJECTS_DIR).
   

#
# Output:
# -------
# - Aligned and processed PET data is stored in the subject's directory within `SUBJECTS_DIR`.
# - PVC and non-PVC results are generated and saved, including cortical surface maps.
# - Intermediate files are cleaned up, leaving only the relevant outputs.
#
# Example:
# --------
#./petsurfer.sh -s subject_label -p pet_dir -t tracer [-r ref_region] [-f fwhm] [-n threads]
#
# This will process the PET data for subject subject_label, using the specified tracer, and perform PVC. 
# The default reference region for SUVR is the cerebellar cortex, 
# the default FWHM for the point spread function is 2.8 and 
# the default number of threads for parallel processing is 1. 
# The results will be saved in the subject's directory within `SUBJECTS_DIR/pet/tracer`.

# Function to display help
show_help() {
  echo "Usage: $(basename $0) -s subject_label -p pet_dir -t tracer [-f fwhm] [-n threads]"
  echo
  echo "Arguments:"
  echo "  -s subject_label: The label of the subject (e.g., sub-01)"
  echo "  -p pet_dir: Path to the PET directory containing the subject data"
  echo "  -t tracer: Tracer name to match in the PET file (e.g., 18Fflortaucipir)"
  echo "  -r ref_region: Reference region for SUVR calculation (e.g., Cerebellum-Cortex)"
  echo "  -f fwhm: (Optional) Full Width at Half Maximum (FWHM) for PSF in PVC"
  echo "  -n threads: (Optional) Number of threads for parallel processing. Default is 1."
  echo "  -h: Show help message"
}

# Default FWHM and threads
FWHM=2.8 # Default FWHM for the point spread function (PSF) in PVC
THREADS=1 # Default number of threads for parallel processing
REF_REGION="Cerebellum-Cortex" # Default reference region for SUVR calculation

# Parse command-line arguments
while getopts "s:p:t:r:f:n:h" opt; do
  case ${opt} in
    s) SUBJECT_LABEL=$OPTARG ;;
    p) PET_DIR=$OPTARG ;;
    t) TRACER=$OPTARG ;;
    r) REF_REGION=$OPTARG ;;
    f) FWHM=$OPTARG ;;
    n) THREADS=$OPTARG ;;
    h) show_help; exit 0 ;;
    *) show_help; exit 1 ;;
  esac
done

# Check required arguments
if [ -z "$SUBJECT_LABEL" ] || [ -z "$PET_DIR" ] || [ -z "$TRACER" ]; then
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
FS_PET_DIR="$SUBJECTS_DIR/$SUBJECT_LABEL/pet/$TRACER"
T1="$SUBJECTS_DIR/$SUBJECT_LABEL/mri/T1.mgz"
TEMPLATE="$FS_PET_DIR/template.mgz"
SURF="$SUBJECTS_DIR/$SUBJECT_LABEL/surf"
PVC="$FS_PET_DIR/pvc"
NOPVC="$FS_PET_DIR/nopvc"
# Path to FreeSurferColorLUT.txt
LUT_FILE="$FREESURFER_HOME/FreeSurferColorLUT.txt"

# Create necessary directories
mkdir -p "$FS_PET_DIR"
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


# Find PET files that match the subject and tracer
PET_FILE=$(find "$PET_DIR/$SUBJECT_LABEL" -type f -name "*trc-${TRACER}_pet.nii.gz" | head -n 1)

if [ -z "$PET_FILE" ]; then
echo "No PET file found for subject $SUBJECT_LABEL with tracer $TRACER in $PET_DIR"
exit 1
fi

echo "Found PET file: $PET_FILE"
# Move and rename the PET file
cp "$PET_FILE" "$FS_PET_DIR/$TRACER.nii.gz"

echo "========Processing PET File==========="
  # Use this PET file for subsequent steps

# Process each PET file
  NUM_FRAMES=$(mri_info "$PET_FILE" --nframes)

  # Align each frame to the T1 volume
  for ((i=0; i<NUM_FRAMES; i++)); do
    FRAME_FILE="$FS_PET_DIR/frame_$i.mgz"
 
    # Extract the i-th frame
    mri_convert "$PET_FILE" \
                --frame $i \
                "$FRAME_FILE"

    # Register the frame to the T1 volume
    mri_coreg --s "$SUBJECT_LABEL" \
              --mov "$FRAME_FILE" \
              --reg "$FS_PET_DIR/frame_$i.reg.lta" \
              --threads $THREADS

    mri_vol2vol --mov "$FRAME_FILE" \
                --targ "$T1" \
                --lta "$FS_PET_DIR/frame_$i.reg.lta" \
                --o "$FS_PET_DIR/frame_reg_$i.mgz"
  done

  # Collect all frame_reg_*.mgz files into a string
  frame_files=$(find "$FS_PET_DIR" -name "frame_reg_*.mgz" | sort | tr '\n' ' ')

  # Ensure that frame_files is not empty before proceeding
  if [ -z "$frame_files" ]; then
    echo "No frame_reg_*.mgz files found in $FS_PET_DIR"
    exit 1
  fi


# Run mri_concat with the collected files
  mri_concat --i $frame_files \
            --o "$FS_PET_DIR/frame_mean.mgz" \
            --mean

# Apply smoothing
  mri_convert --fwhm $FWHM \
              "$FS_PET_DIR/frame_mean.mgz" \
              "$TEMPLATE" \
              --nthreads $THREADS

# Obtain the brain mask
  mri_mask "$TEMPLATE" \
          "$SUBJECTS_DIR/$SUBJECT_LABEL/mri/brainmask.mgz" \
          "$FS_PET_DIR/brainmasked.mgz" || exit 1

echo "========Coregistration==========="
# Register the template to generate a template.reg.lta file
  mri_coreg --s "$SUBJECT_LABEL" \
            --mov "$TEMPLATE" \
            --reg "$FS_PET_DIR/template.reg.lta" \
            --threads $THREADS

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
  
echo "========Partial Volume Correction (PVC)==========="

  # Perform Partial Volume Correction (PVC)
  mri_gtmpvc --i "$TEMPLATE" \
            --reg "$FS_PET_DIR/template.reg.lta" \
            --psf $FWHM \
            --seg gtmseg.mgz \
            --default-seg-merge \
            --auto-mask 1 .05 \
            --mgx .25 \
            --rescale $REF_REGION_VALUES \
            --save-input \
            --o "$PVC" \
            --threads $THREADS || exit 1

echo "========Mapping PVC to Surface==========="

  mri_vol2surf --src "$PVC/mgx.ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi lh \
              --projfrac 1 \
              --o  "$SURF/lh_pet_pvc_$TRACER.nii.gz" \
              --no-reshape || exit 1

  mri_vol2surf --src "$PVC/mgx.ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi rh \
              --projfrac 1 \
              --o  "$SURF/rh_pet_pvc_$TRACER.nii.gz" \
              --no-reshape || exit 1

echo "==Calculating SUVR without Partial Volume Correction (noPVC)=="

  # Without Partial Volume Correction (noPVC)

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
              --o  "$SURF/lh_pet_nopvc_$TRACER.nii.gz" \
              --no-reshape || exit 1

  mri_vol2surf --src "$NOPVC/ctxgm.nii.gz" \
              --srcreg "$PVC/aux/bbpet2anat.lta" \
              --hemi rh \
              --projfrac 1 \
              --o  "$SURF/rh_pet_nopvc_$TRACER.nii.gz" \
              --no-reshape


echo "========Converting Surfaces to Curv==========="

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/lh_pet_pvc_$TRACER.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $SURF/pet_pvc_$TRACER \
  --hemi lh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/rh_pet_pvc_$TRACER.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $SURF/pet_pvc_$TRACER \
  --hemi rh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/lh_pet_nopvc_$TRACER.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $SURF/pet_nopvc_$TRACER \
  --hemi lh \
  --trg_type curv

  mri_surf2surf --srcsubject $SUBJECT_LABEL \
  --srcsurfval $SURF/rh_pet_nopvc_$TRACER.nii.gz \
  --trgsubject $SUBJECT_LABEL \
  --trgsurfval $SURF/pet_nopvc_$TRACER \
  --hemi rh \
  --trg_type curv

echo "========Auxiliary File Cleanup==========="
# Clean up frame files
# Collect all frame_reg_*.mgz files into a string

frame_files=$(find "$FS_PET_DIR" -name "frame_reg_*.mgz" | sort | tr '\n' ' ')
  rm $frame_files
find "$FS_PET_DIR" -name "frame_*" ! -name "frame_mean.mgz" -exec rm -f {} \;

echo "========Processing completed for subject $SUBJECT_LABEL==========="
echo "Surface projections: $SURF/pet_nopvc_$TRACER and $SURF/pet_pvc_$TRACER"
