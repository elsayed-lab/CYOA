package Bio::Adventure::Proteomics;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

=head1 NAME

Bio::Adventure::Proteomics - Invoke tools to process proteomics data, mostly DIA-SWATH.

=head1 SYNOPSIS

=head1 METHODS

=cut

## Since I have not been doing any of these analyses in a while, I want to start by
## just getting my invocations written down in this file.  This is mostly copy/pasted
## from my interactive shell script 'dia_invocation_20191001.sh'.
## It in turn sets quite a few options from the file 'parameters/20191001_settings.sh'.

sub Parameters {

    my $jstring = qq!
export START="/fs/cbcb-scratch/abelew/mycobacterium_tuberculosis_2019"

echo "Setting variables for the rest of the script."
export MZ_WINDOWS="${MZ_WINDOWS:-8mz}"
export VERSION="${VERSION:-201805}"
export TYPE="${TYPE:-whole}"
export SEARCH_METHOD="${SEARCH_METHOD:-comet}"
export DDA_METHOD="${DDA_METHOD:-HCD}"
export TEST="${TEST:-Rv1611}"
export CLEAVAGES="${CLEAVAGES:-2}"
export REFDB="${REFDB:-reference/reference.fasta}"
export DECOY_STRING="${DECOY_STRING:-DECOY_}"
export CONFIDENCE_CUTOFF="${CONFIDENCE_CUTOFF:-0.05}"
export IPPP="${IPPP:-0.95}"

## Parameters from here on depend on previous parameters
export MZXML_DIR="preprocessing/01mzXML"
##export DDA_INPUTS=$(/bin/ls "${MZXML_DIR}/dda/${VERSION}/"*.mzXML)
##export PEPXML_INPUTS=$(/bin/ls "${MZXML_DIR}/dda/${VERSION}/"*.pep.xml)

export FDR_DIR="preprocessing/02fdr_controlled/${VERSION}"
mkdir -p "${FDR_DIR}"
export FDR_RESULT="${FDR_DIR}/fdr_library_${TYPE}_${SEARCH_METHOD}_${DDA_METHOD}.xml"

export IPP_DIR="preprocessing/03iprophet/${VERSION}"
mkdir -p "${IPP_DIR}"
export IPP_RESULT="${IPP_DIR}/iprophet_${TYPE}_${SEARCH_METHOD}_${DDA_METHOD}.xml"

export MAYU_DIR="preprocessing/04mayu/${VERSION}"
mkdir -p "${MAYU_DIR}"
export MAYU_RESULT="${MAYU_DIR}/mayu_${TYPE}_${SEARCH_METHOD}_${DDA_METHOD}.out"

export SPEC_DIR="preprocessing/05spectral_libraries/${VERSION}"
mkdir -p "${SPEC_DIR}"
export SPECTRAL_BASENAME="${SPEC_DIR}/${TYPE}_${SEARCH_METHOD}_${DDA_METHOD}"
export SPECTRAL_LIBRARY="${SPECTRAL_BASENAME}.splib"

export CON_DIR="preprocessing/06spectral_consensus/${VERSION}"
mkdir -p "${CON_DIR}"
export CONSENSUS="${CON_DIR}/${TYPE}_${SEARCH_METHOD}_${DDA_METHOD}"
export CONSENSUS_TRAML="${CONSENSUS}_${MZ_WINDOWS}.TraML"
export OPTIMIZED_TRAML="${CONSENSUS}_${MZ_WINDOWS}_optimized.TraML"
export DECOY_TRAML="${CONSENSUS}_${MZ_WINDOWS}_decoy.TraML"
export TRANSITION_PREFIX="${CONSENSUS}_${MZ_WINDOWS}_transitions"
export TUBERCULIST_PQP="preprocessing/spectral_consensus/tuberculist/201805_tuberculist_comet_HCD_20mz_decoy.pqp"

export OSW_DIR="preprocessing/07openswath/${VERSION}"
mkdir -p "${OSW_DIR}"
export SWATH_OUTDIR="${OSW_DIR}/${TYPE}_${MZ_WINDOWS}"
export TUBERCULIST_OUTDIR="${SWATH_OUTDIR}_tuberculist"

export PYP_DIR="preprocessing/08pyprophet/${VERSION}"
mkdir -p "${PYP_DIR}"
export PYPROPHET_OUTDIR="${PYP_DIR}/${TYPE}_${MZ_WINDOWS}"

export TRIC_DIR="preprocessing/09tric/${VERSION}"
export TRIC_OUTDIR="${TRIC_DIR}/${TYPE}_${MZ_WINDOWS}"
mkdir -p "${TRIC_OUTDIR}"

echo "Loading environment modules."
module add openms

!;

}

## Invoke the DDA search tool, comet and use refreshparser in order to add back
## potential alternate matches which were filtered out by comet.
sub Comet {

    my $jstring = qq!comet
  -Pparameters/comet_${DDA_METHOD}_params.txt \\
  ${comet_input} \\
  2>>"${comet_input}.log" 1>&2
RefreshParser \\
  ${refresh_input} \\
  "${REFDB}" \\
  2>"${refresh_input}_refreshparser.out" 1>&2

!;

}

## Use xinteract to merge multiple DDA searches
sub Merge_Xinteract {

    my $jstring = qq!
xinteract \\
    "-d${DECOY_STRING}" \\
    -OARPpd \\
    -Ow \\
    "-N${FDR_RESULT}" \\
    ${PEPXML_INPUTS} \\
    2>"${FDR_RESULT}.log" 1>&2
!;

}

## Use interprophetparser to combine identifications.

sub Merge_Peptide_Identifications {

    my $jstring = qq!
InterProphetParser \\
    "DECOY=${DECOY_STRING}" \\
    "${FDR_RESULT}" \\
    "${IPP_RESULT}" \\
    2>"${IPP_RESULT}.log" 1>&2
!;
}

sub Standardize_FDR_Mayu {

    my $jstring = qq!
Mayu.pl \\
    -A "${IPP_RESULT}" \\
    -C "${REFDB}" \\
    -E "${DECOY_STRING}" \\
    -G "${CONFIDENCE_CUTOFF}" \\
    -H 101 \\
    -I "${CLEAVAGES}" \\
    2>"${MAYU_RESULT}" 1>&2
!;
}

sub Create_Spectral_Libraries {

    my $jstring = qq?
spectrast \\
    "-cN${SPECTRAL_BASENAME}" \\
    "-cI${DDA_METHOD}" \\
    -cf "Protein! ~ ${DECOY_STRING}" \\
    "-cP${IPPP}" \\
    -c_IRR \\
    "${IPP_RESULT}" \\
    2>"${IPP_RESULT}.out" 1>&2
?;
}

sub Make_Consensus_Libraries {

    my $jstring = qq?
spectrast \\
    "-cN${CONSENSUS}" \\
    "-cI${DDA_METHOD}" \\
    -cAC \\
    "${SPECTRAL_LIBRARY}" \\
    2>"${SPECTRAL_LIBRARY}.out" 1>&2
?;
}

sub Create_Acquisition_Windows {

    my $jstring = qq!
echo "Make sure you have created an acquisition window file _without_ a header."
echo "There should be some which are just bob.txt in the windows/ directory."
echo "Writing spectrast tsv files with spectrast2tsv."
echo "Making a consensus library specific for ${MZ_WINDOWS} windows."
echo "The input is: ${CONSENSUS}.sptxt"
echo "The output is: ${CONSENSUS}_${MZ_WINDOWS}.tsv"
spectrast2tsv.py \\
    -l 350,2000 \\
    -s b,y \\
    -x 1,2 \\
    -o 4 \\
    -n 6 \\
    -p "${CONFIDENCE_CUTOFF}" \\
    -d \\
    -e \\
    -w "windows/acquisition_${MZ_WINDOWS}.txt" \\
    -k openswath \\
    -a "${CONSENSUS}_${MZ_WINDOWS}.tsv" \\
    "${CONSENSUS}.sptxt" \\
    2>"${CONSENSUS}_spectrast.log" 1>&2
!;
}

sub Spectral_Convert {

    my $jstring = qq!
echo "Converting spectral libraries to TraML."
rm -f "${CONSENSUS}_${MZ_WINDOWS}.TraML"
echo "The input is: ${CONSENSUS}_${MZ_WINDOWS}.tsv"
echo "The output is: ${CONSENSUS}_${MZ_WINDOWS}.TraML"
TargetedFileConverter \\
    -in "${CONSENSUS}_${MZ_WINDOWS}.tsv" \\
    -in_type tsv \\
    -out "${CONSENSUS}_${MZ_WINDOWS}.TraML" \\
    -out_type TraML

echo "Converting spectral libraries to TraML."
rm -f "${CONSENSUS}_${MZ_WINDOWS}.TraML"
echo "The input is: ${CONSENSUS}_${MZ_WINDOWS}.tsv"
echo "The output is: ${CONSENSUS}_${MZ_WINDOWS}.TraML"
TargetedFileConverter \\
    -in "${CONSENSUS}_${MZ_WINDOWS}.tsv" \\
    -in_type tsv \\
    -out "${CONSENSUS}_${MZ_WINDOWS}.TraML" \\
    -out_type TraML

echo "Converting the Tuberculist libraries to pqp."
TUBERCULIST_TRAML="results/05spectral_libraries/Mtb_TubercuList-R27_iRT_UPS_noMox_noMC_sall_osw_decoy.TraML"
echo "The input is: ${TUBERCULIST_TRAML}"
echo "The output is: ${TUBERCULIST_PQP}"
TargetedFileConverter \\
    -in "${TUBERCULIST_TRAML}" \\
    -in_type TraML \\
    -out "${TUBERCULIST_PQP}" \\
    -out_type pqp \\
    2>"${TUBERCULIST_PQP}_convert.log" 1>&2

!;

}

sub Optimize_Spectral_Libraries {

    my $jstring = qq!
## Apparently the append/exclude_similar/etc options were removed from the newer
## version of openms.
echo "Adding decoys and optimizing the spectral libraries."
rm -f "${CONSENSUS}_${MZ_WINDOWS}_decoy.TraML"
echo "Generating decoys for ${CONSENSUS}_${MZ_WINDOWS}"
##-min_transitions 6 \
##-max_transitions 6 \
##-allowed_fragment_types b,y \
##-allowed_fragment_charges 1,2,3,4 \
##-enable_detection_specific_losses \
##-enable_detection_unspecific_losses \
##-precursor_mz_threshold 0.025 \
##-precursor_lower_mz_limit 400 \
##-precursor_upper_mz_limit 1200 \
##-product_mz_threshold 0.025 \
##-product_lower_mz_limit 350 \
##-product_upper_mz_limit 2000 \
echo "Optimizing the consensus library."
echo "The input is: ${CONSENSUS_TRAML}"
echo "The output is: ${OPTIMIZED_TRAML}"
OpenSwathAssayGenerator \\
     -in "${CONSENSUS_TRAML}" \\
     -out "${OPTIMIZED_TRAML}" \\
     -swath_windows_file "windows/openswath_${MZ_WINDOWS}.txt" \\
     -enable_ipf \\
     -unimod_file "parameters/unimod.xml" \\
     2>"${OPTIMIZED_TRAML}.log" 1>&2
echo "Adding decoys."
echo "The input is: ${OPTIMIZED_TRAML}"
echo "The output is: ${DECOY_TRAML}"
OpenSwathDecoyGenerator \\
    -in "${OPTIMIZED_TRAML}" \\
    -out "${DECOY_TRAML}" \\
    -decoy_tag "${DECOY_STRING}" \\
    -method shuffle \\
    2>"${DECOY_TRAML}.log" 1>&2
grep "${TEST}" "${DECOY_TRAML}" | head

echo "Converting decoy-added libraries back to tsv for examination later."
rm -f "${CONSENSUS}_${MZ_WINDOWS}_decoy.tsv"
echo "The input is: ${DECOY_TRAML}"
echo "The output is: ${TRANSITION_PREFIX}.tsv"
TargetedFileConverter \\
    -in "${DECOY_TRAML}" \\
    -in_type TraML \\
    -out "${TRANSITION_PREFIX}.tsv" \\
    -out_type tsv \\
    2>"${TRANSITION_PREFIX}_convert_tsv.log" 1>&2
echo "Converting libraries with decoys for ${TRANSITION_PREFIX} to pqp."
echo "The input is: ${DECOY_TRAML}"
echo "The output is: ${TRANSITION_PREFIX}.pqp"
TargetedFileConverter \\
    -in "${TRANSITION_PREFIX}.tsv" \\
    -in_type tsv \\
    -out_type pqp \\
    -out "${TRANSITION_PREFIX}.pqp"
echo "Transition library for openswathworkflow is: ${TRANSITION_PREFIX}.pqp"
grep "${TEST}" "${TRANSITION_PREFIX}.tsv" | head
grep -c "${DECOY_STRING}" "${TRANSITION_PREFIX}.tsv"

!;
}

sub Openswath_vs_Comet_Transitions {

    my $jstring = qq!
## Here is a partial explanation of some options:
## https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=5cda21e7a2eee96f&render_repository_actions_for=tool_shed&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F002%2Frepo_2999%2FOpenSwathWorkflow.xml&changeset_revision=380ced982ec7
## -rt_extraction_window: Only extract RT around this value (-1 means extract over the whole range, a value of 600 means extract around +/0 300 of the expected elution.)
## -mz_extraction_window: Size of window, add the -ppm flag to use parts per million.
## -rt_normalization_factor: The normalized RT is expected to go from 0-1, if it goes outside this range, modidify this option.
## -uis_threshold_sn: S/N threshold to consider identification transition (-1 looks at all)
## -min_peak_width: Discard peaks with a width < than this, -1 discards no peaks.
## -background_subtraction: background is normally estimated at peak boundaries, this can change that behavior
## -recalculate_peaks: Tries to use the consensus (median) peak border if the variation of picked peaks is too large.
## -recalculate_peaks_max_z: If the peak is beyond this z-score, then recalculate using the median as above.
## -minimal_quality: if compute_peak_quality is set, then this sets the lower boundary of the quality score.
## -compute_peak_quality: Calculate a quality score centered on 0 where good is >= 0, bad is <= -1.
## -sgolay_frame_length: Number of data points used for smoothing, it must be an odd number.
## -gauss_width: Estimated peak width in gauss units.
## -use_gauss: Use a gaussian filter for smoothing instead of the default savitzky-golay filter.
## -peak_width: Force a minimum peak width on the data by extending by this amount, -1 turns this off.
## -signal_to_noise: threshold at which peaks will not be extended further, too high is problematic.
## -write_sn_log_messages: Write out signal noise messages, I wish it wouldn't.
## -remove_overlapping_peaks: Remove overlapping peaks during picking (ooohh I should turn this on)
## -method: Method to choose for chromatographic peak-picking.
## -dia_extraction_window: By default Th, -ppm will switch it.
## -dia_centroided: Use centroided data.
## -dia_byseries_intensity_min: b/y minimum intensities to consider.
## -dia_byseries_ppm_diff: Minimum b/y difference in ppm to consider.
## -dia_nr_isotopes: NR of isotopes to consider.
## -peak_before_mono_max_pp_diff: Maximum ppm difference to count a peak at lower m/z while looking for non-monoisotopic peaks.
## -max_iteration: Maximum number of iterations when performing the Levenberg-Marquardt algorithm.

echo "Invoking the OpenSwathWorkflow using local comet-derived transitions."
echo "Checking in, the transition library is: ${TRANSITION_PREFIX}.pqp"
swath_inputs=$(/bin/ls "results/mzXML/dia/${VERSION}/")
echo "Checking in, the inputs are: ${swath_inputs}"
mkdir -p "${SWATH_OUTDIR}"
mkdir -p "${SWATH_OUTDIR}_tuberculist"
for input in ${swath_inputs}
do
    type="comet"
    in_mzxml="results/mzXML/dia/${VERSION}/${input}"
    name=$(basename "${input}" .mzXML)
    echo "Starting ${type} libraries run of ${name} using ${MZ_WINDOWS} windows at $(date)."
    swath_output_prefix="${SWATH_OUTDIR}_${type}/${name}_${DDA_METHOD}"
    pyprophet_output_prefix="${PYPROPHET_OUTDIR}_${type}/${name}_${DDA_METHOD}"
    echo "Deleting previous swath output file: ${swath_output_prefix}.osw"
    rm -f "${swath_output_prefix}.osw"
    OpenSwathWorkflow \\
        -ini "parameters/openms_${VERSION}.ini" \\
        -in "${in_mzxml}" \\
        -swath_windows_file "windows/openswath_${name}.txt" \\
        -tr "${TRANSITION_PREFIX}.pqp" \\
        -out_osw "${swath_output_prefix}.osw" \\
        2>"${swath_output_prefix}_osw.log" 1>&2
    echo "Scoring individual swath run: ${swath_output_prefix}"
    pyprophet \\
        score \\
        --level ms1 \\
        --in "${swath_output_prefix}.osw" \\
        --out "${pyprophet_output_prefix}.osw" \\
        2>>"${pyprophet_output_prefix}_ms1.log" 1>&2
    pyprophet \\
        score \\
        --level ms2 \\
        --in "${pyprophet_output_prefix}.osw" \\
        --out "${pyprophet_output_prefix}.osw" \\
        2>>"${pyprophet_output_prefix}_ms2.log" 1>&2
    pyprophet \\
        protein \\
        --in "${pyprophet_output_prefix}.osw" \\
        --context run-specific \\
        2>>"${pyprophet_output_prefix}_protein.log" 1>&2
    echo "Exporting individual swath run: to ${pyprophet_output_prefix}.tsv"
    pyprophet \\
        export \\
        --in "${pyprophet_output_prefix}.osw" \\
        --out "${pyprophet_output_prefix}.tsv" \\
        2>>"${pyprophet_output_prefix}_export.log" 1>&2
done

type="our"
echo "Creating final aligned matrix for the ${type} library."
echo "The inputs are the tsv files in: ${SWATH_OUTDIR}_${type}"
echo "The output is: ${TRIC_OUTDIR}_${lib}/${SEARCH_METHOD}_${DDA_METHOD}.tsv"
feature_alignment.py \\
    --force \\
    --target_fdr "${CONFIDENCE_CUTOFF}" \\
    --max_fdr_quality "${CONFIDENCE_CUTOFF}" \\
    --in ./${SWATH_OUTDIR}_${lib}/*.tsv \\
    --out "${TRIC_OUTDIR}_${lib}/${SEARCH_METHOD}_${DDA_METHOD}.tsv" \\
    --out_matrix "${TRIC_OUTDIR}_${lib}/${DDA_METHOD}_outmatrix.tsv" \\
    --out_meta "${TRIC_OUTDIR}_${lib}/${DDA_METHOD}_meta.tsv" \\
    2>"${TRIC_OUTDIR}_${lib}/${SEARCH_METHOD}_${DDA_METHOD}.log" 1>&2
tail "${TRIC_OUTDIR}_${lib}/${SEARCH_METHOD}_${DDA_METHOD}.log"

!;

}

sub Openswath_vs_Sequence_Transitions {

    my $jstring = qq!
echo "Invoking the OpenSwathWorkflow using the tuberculist transitions."
base_mzxmldir="results/01mzXML/dia/${VERSION}"
swath_inputs=$(/bin/ls "${base_mzxmldir}")
echo "Checking in, the inputs are: ${swath_inputs}"
mkdir -p "${TUBERCULIST_OUTDIR}"
pypdir="${PYPROPHET_OUTDIR}_tuberculist"
mkdir -p "${pypdir}"
for input in ${swath_inputs}
do
    in_mzxml="${base_mzxmldir}/${input}"
    name=$(basename "${input}" .mzXML)
    echo "Starting openswath run of ${name} using ${MZ_WINDOWS} windows at $(date)."
    tb_output_prefix="${TUBERCULIST_OUTDIR}/${name}_vs_${VERSION}_${TYPE}_${DDA_METHOD}_dia"
    pyprophet_output_prefix="${pypdir}/${name}_vs_${VERSION}_${TYPE}_${DDA_METHOD}_dia"
    echo "Deleting previous swath output file: ${tb_output_prefix}.osw"
    rm -f "${tb_output_prefix}.osw"
    OpenSwathWorkflow \\
        -ini "parameters/openms_${VERSION}.ini" \\
        -in "${in_mzxml}" \\
        -swath_windows_file "windows/openswath_${name}.txt" \\
        -tr "${TUBERCULIST_PQP}" \\
        -out_osw "${tb_output_prefix}.osw" \\
        2>"${tb_output_prefix}_osw.log" 1>&2
    if [[ "$?" -ne "0" ]]; then
        echo "OpenSwathWorkflow for ${name} failed."
    fi

    rm -f "${tb_output_prefix}_scored.osw"
    echo "Scoring individual swath run: ${tb_output_prefix}"
    pyprophet \\
        score \\
        --level ms1 \\
        --in "${tb_output_prefix}.osw" \\
        --out "${pyprophet_output_prefix}_scored.osw" \\
        2>>"${pyprophet_output_prefix}_pyprophet_ms1.log" 1>&2
    if [[ "$?" -ne "0" ]]; then
        echo "MS1 scoring ${pyprophet_output_prefix}_scored.osw failed."
    fi

    pyprophet \\
        score \\
        --level ms2 \\
        --in "${pyprophet_output_prefix}_scored.osw" \\
        --out "${pyprophet_output_prefix}_scored.osw" \\
        2>>"${pyprophet_output_prefix}_pyprophet_ms2.log" 1>&2
    if [[ "$?" -ne "0" ]]; then
        echo "MS2 scoring ${pyprophet_output_prefix}_scored.osw failed."
    fi

    pyprophet \\
        protein \\
        --in "${pyprophet_output_prefix}_scored.osw" \\
        --context run-specific \\
        2>>"${pyprophet_output_prefix}_pyprophet_protein.log" 1>&2
    if [[ "$?" -ne "0" ]]; then
        echo "Protein scoring ${pyprophet_output_prefix}_scored.osw failed."
    fi

    rm -f "${pyprophet_output_prefix}_scored.tsv"
    echo "Exporting individual swath run: to ${pyprophet_output_prefix}_scored.tsv"
    pyprophet \\
        export \\
        --in "${pyprophet_output_prefix}_scored.osw" \\
        --out "${pyprophet_output_prefix}_scored.tsv" \\
        2>>"${pyprophet_output_prefix}_pyprophet_export.log" 1>&2
    ## ok something is fubar, the stupid tsv files are being written in the cwd as run_filename.tsv
    ## No matter what I do
    mv "${input}.tsv" "${pyprophet_output_prefix}_scored.tsv"
    if [[ "$?" -ne "0" ]]; then
        echo "Exporting ${pyprophet_output_prefix}_scored.tsv failed."
    fi
done

tric_tb="${TRIC_OUTDIR}_tuberculist"
mkdir -p "${tric_tb}"
feature_alignment.py \\
    --force \\
    --in "./${pypdir}/"*.tsv \\
    --out "${tric_tb}/${SEARCH_METHOD}_${DDA_METHOD}.tsv" \\
    --out_matrix "${tric_tb}/${DDA_METHOD}_outmatrix.tsv" \\
    --out_meta "${tric_tb}/${DDA_METHOD}_meta.tsv" \\
 2>"${tric_tb}/feature_alignment.err" \\
 1>"${tric_tb}/feature_alignment.out"
echo "Wrote final output to ${tric_tb}/${SEARCH_METHOD}_${DDA_METHOD}.tsv"

!;

}

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure>

=cut

1;
