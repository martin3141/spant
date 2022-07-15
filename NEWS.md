# spant 2.1.0
* y-axis label can now be changed with a new argument (yaxis_lab) to
plot.mrs_data.
* Segmentation files are now checked for non-integer values.
* gridplot has better row and col number defaults for MRSI, and has a bug fix
for plotting masked spectra.
* Updated functions for dealing with lists of MRSI datasets.

# spant 2.0.0
* magrittr pipe (%>%) is no longer loaded together with spant now the base pipe
is part of R. Vignettes and tests have been updated to use the base pipe.
* NIfTI MRS writer no longer fails when dealing with a NULL affine matrix.
* Bug fix for NIfTI MRS writer, where frequency domain data was not correctly
converted to time-domain before export.

# spant 1.22.0
* Added MSM and acetone to the set of internally defined molecules.
* write_mrs_nifti function removed from userspace. write_mrs should be used
instead.
* set_td_pts and crop_td_pts now work with lists of spectra.
* Multiple files can now be written by passing a list of spectra to write_mrs.

# spant 1.21.0
* Bug fix for get_mrsi_voi function due to incorrect rounding of MRSI matrix
dimensions in some cases.
* Added mask_fit_res function.
* Option added to get_2d_psf function to ensure odd output dimensions.
* Bug fix for using rats on datasets with masked voxels.
* Added spec_decomp function.
* Added ret_scale_factor option to scale_spec function.

# spant 1.20.0
* Options added to ft_dyns to return the modulus of the 2D spectra and to
pre-transform the chemical shift axis if needed.
* Bug fix for printing lists of mrs_data objects.
* comb_coils function now supports lists of mrs_data objects.
* IMA reader now respects the verbose option in read_mrs.
* Bug fix for get_mrsi_voi ignoring the ker argument.
* Bug fix for reading Siemens TWIX files.

# spant 1.19.0
* Added sum_mrs_list function.
* Bug fix for operating on mrs_data with single values that aren't numeric
types. For example, multiplying mrs_data by a single complex number now works.
* A number of common preprocessing functions now operate on lists of mrs_data
objects.
* Change to NIfTI MRS reader to fix a failing test when using the development
version of RNifti.

# spant 1.18.0
* Better behaviour when using stackplot with single spectra datasets.
* Better behaviour when using image with single spectra datasets.
* Bug fix for reading MEGA-PRESS data in TWIX format.
* Better error checking for lists of spectra.
* Bug fix for calculating NAA linewidth from MEGA-PRESS data.

# spant 1.17.0
* New function added (recon_twix_2d) for reconstructing basic phase encoded twix
2D MRSI data.
* File name argument for read_mrs now permits globbing, eg read_mrs("*.dcm").
* Improved plotting of metabolite maps containing infinite values.
* Improvements to GE p-file reader.
* Bug fix for TWIX MRSI voxel dimensions.
* get_mrsi2d_seg now returns partial volume maps as well as a data frame.
* Improved error handling for spec_op function.
* apodise_xy can now be applied to multi-coil and dynamic datasets.

# spant 1.16.0
* The package has been published in the Journal for Open Source Software :
"spant: An R package for magnetic resonance spectroscopy analysis. JOSS 2021,
6(67), 3646".
* The basis result of the HSVD function is now an mrs_data object.
* The complete model signal is now output by the HSVD function.
* image function x-axis updated to be consistent with other plotting methods.
* Minor refactor of the simulation code and a dependency swap from complexplus
to exmp packages.

# spant 1.15.0
* Added a unit test for reading and writing LCM .basis formatted files.
* Added FWHM estimates for tCr and tCho in ABfit.
* Added water suppression efficiency and water FWHM measures to the diagnostic
output of svs_1h_brain_analysis.
* Added get_head_dyns and get_tail_dyns to return the first and last dynamic
scans within a dataset.
* Fixed CI errors.
* Improved installation instructions.
* Removed comb_csv_results function and reduced the number of required packages.

# spant 1.14.0
* Added glycerol simulation parameters, e.g. get_mol_paras("glyc").
* Bug fix for read_ima_* functions.
* Improved y = 0 baseline for stackplot when setting bl_lty parameter.
* Removed norm_mrs function and replaced with scale_spec for simple data
scaling tasks, eg scaling based on the integration of a spectral region.
* Added spec_op function for performing simple summary operations on spectral
segments.
* Changed the name of scale_mrs function to scale_mrs_amp.
* Added mean_mrs_list function.
* Improved LCM RAW and BASIS readers (contribution from Alex Craven).

# spant 1.13.0
* ABfit frequency shifts limits are now specified in ppm rather than Hz to
improve consistency between field strengths.
* NAA linewidth is now estimated and output by ABfit when NAAG is absent - 
useful for BRAINO phantom scans.
* Warning now given when spectra are mathematically combined and are not both
in the same time/frequency domain.
* sum_mrs function added to simplify combining spectra within a pipe.
* scale_mrs function added to simplify scaling a spectrum within a pipe.
* add_noise function added to simplify generating simulated data with a pipe.
* hsvd_filt function now accepts a frequency range in ppm units and gives the
option to return the model rather than the filtered data.
* Bug fix for hsvd_filt function where the max_damp argument was ignored.
* SVS reference scans, found in some TWIX files, are now removed by default.
* Bug fix for reading list/data when only one coil element is used.

# spant 1.12.0
* Added 2HG and citrate simulation parameters. e.g. get_mol_paras("2hg").
* Better print output for molecular definitions.
* Added metabolite and basis simulation vignettes.
* Bug fix for setting the ppm reference when reading LCModel RAW files.
* Bug fix for ABfit CRLB calculation of combined signals, eg tNAA, tCr.

# spant 1.11.0
* Options added to allow extra information to attached to mrs_data and
fit_result objects as a data frame.
* New functions (combine_fit_XXX) for working with multiple fit results
contained within a list structure.
* Basis set and noise region are now checked for validity in ABfit.
* Added spant_abfit_benchmark function.
* Tentative functions for performing "standard" 1H brain analyses: 
svs_1h_brain_analysis and svs_1h_brain_batch_analysis.
* Improved support for LCModel analyses.

# spant 1.10.0
* Fix for NIfTI MRS reader/writer.
* ortho3 now shows correct labels for orientations other than RAS.
* ortho3_int function renamed to ortho3_inter.
* Argument order change to plot_voi_overlay and plot_voi_overlay_seg to be more
consistent with ortho3.
* Regression fix for partial volume segmentation plotting.
* Echo time parameter is now stored in the meta structure.
* ABfit now performs a 1D phase parameter search before the prefit stage to
improve reliability. May be disabled with the prefit_phase_search fit option.

# spant 1.9.0
* NIfTI MRS reader and writer now uses the header extension for metadata. Thanks
to Jon Clayden for adding extension read/write support to the RNifti package.
* Default plots now have gridlines in the y-direction and the plot line is now
thicker and coloured blue.
* Opacity option added to the plotting functions (alpha).
* Bug fix for comb_coils with SVS data.
* Bug fixes for Siemens geometry information.
* Changed the ordering of arguments to write_mrs and write_mrs_nifti to improve 
consistency with other functions.
* Internal function (ortho3) now used to plot voxel locations on MRI.
* Internal dicom reader function added.
* Tentative support for Siemens and Philips DICOM MRS format added.
* Internal changes to the way orientation information is handled.

# spant 1.8.0
* Added gridplot function.
* New functions added for down sampling.
* Added signal space projection method for MRSI.
* Geometry information is now read from Siemens twix files.
* GitHub actions are now used for continuous integration instead of Travis and
AppVeyor.
* Added precomp function to avoid repeated computation.
* mrs_data2mat function now collapses all dimensions to dynamics by default.
* mrs_data objects now store the nucleus.
* Added a reader for old Varian format (fid/procpar) data.

# spant 1.7.0
* Added write_mrs function which guesses the output format from the file
extension or can be specified as an argument. write_mrs_XXX functions have been
depreciated.
* read_mrs function now tries to guess the format from the file extension.
* Added json sidecar to NIFTI MRS export function.
* Added the option to read MRS data from a NIFTI file and json sidecar using the
read_mrs function.
* Changed default crop_spec region to between 4.0 and 0.2 ppm.

# spant 1.6.0
* Bug fix for GE P file reader.
* Added downsample_mrs function.
* Tentative function (write_mrs_nifti) to write MRS data as a NIFTI format file 
- for evaluation purposes only.

# spant 1.5.0
* Added an option to ABfit to allow the metabolite amplitudes to be negative 
(ahat_calc_method).
* Removed lsei package dependency.

# spant 1.4.0
* Added the option to plot a y = 0 baseline trace for stackplot.mrs_data.
* Added convenience functions to read and write mrs_data to rds format.

# spant 1.3.0
* Added get_fit_table function to combine all fit tables in a fit_result object
into a single dataframe.
* ppm function can now be applied to fit result objects.
* Bug fix for plot_slice_fit_inter and added the option to specify a
denominator.
* Added the option to specify a denominator to plot_slice_fit.

# spant 1.2.1
* Added the option to display a progress bar in fit_mrs function for better
conformance to "Writing R Extensions" in non-interactive use.
* Changed test tolerance to accommodate differences with OSX.

# spant 1.2.0
* Performance improvement for HSVD water filter.
* ABfit unit tests are now run on simulated data to improve consistency between
different platforms.
* Improvements to fit amplitude scaling code.
* Improved checking for mrs_data processing functions.
* Added preprocessing steps vignette.

# spant 1.1.0
* Improved ppm labels for ABfit plot results.
* Bug fix for plot_slice_fit when using fits from masked data.
* Updated unit tests.
* Added vignette on manually adjusting ABfit baseline smoothness.
* Updated the introduction vignette.
* Bug fixes for image and stackplot functions for masked MRS data.
* Bug fix for using RATS with masked MRS data.
* Changes to prepare for for R 4.0.0.

# spant 1.0.0
* ABfit analysis method has been added, and is now default for the mrs_fit
function.
* Added reader for LCModel RAW format data.
* Added read_mri_dyn_dir function for reading dynamic MRS exported from Siemens
scanners.
* Bug fixes for 2D MRSI voxel segmentation calculation.
* Bug fixes for Siemens IMA format reader for SVS data.
* Optional colourbar added to ortho3 function.

# spant 0.19.0
* Added Asc, BHB, Cho, PEth and Ser simulation parameters.
* Added ker option to get_mrsi_voi function.
* Added append_basis function to combine two basis sets objects.
* The align function now accepts more than one reference frequency.

# spant 0.17.0
* Added a function to grid shift 2D MRSI data in the x/y direction.
* Better plotting/fitting support for masking data by setting data points to NA.
* Bug fix for interactive voxel selection position indicator.
* Added mask_xy to mask voxels in a centred rectangular region.
* Minor changes to improve parallel processing support.

# spant 0.16.0
* SNR is now defined as the max data point divided by the standard deviation of
the noise (n.b. factor of two has been removed in-line with upcoming terminology
paper).
* Default rats method improved to work with multidimensional datasets.
* Added norm_mrs function to normalise the intensity of spectral data.
* Added bc_constant function to correct spectral baselines by a constant offset.
* Added re_weighting function to apply a resolution enhancement weighting to the
FID.
* Performance improvement for apodise_xy function.
* sd function now works for mrs_data.
* Added 2D MRSI support for Siemens IMA format.

# spant 0.15.0
* Bug fix for using auto_phase function with a single spectrum.
* Bug fix for comb_coils not returning unaveraged data when requested.
* Added options to combine metabolite signals from the stackplot of a fit object
and adjust the plot margins.
* Added comb_fits function.
* Added collapse_to_dyns function for mrs_fit objects.
* RDA reader now extracts geometry information.

# spant 0.14.0
* Added options to omit basis signals, change label names and combine lipid and
MM signals from the stackplot of a fit object.
* Added auto_phase function for zeroth order phase-correction of simple spectra.
* Added get_subset function to aid MRSI slicing.
* Added decimate_mrs function.
* Added fit_amps function to quickly extract amplitude estimates from a fit
object.
* Bug fix for int_spec function.
* sim_basis function arguments updated to accept acq_par objects.

# spant 0.13.0
* Various bug fixes for Siemens TWIX reader.
* rats and tdsr functions now use the mean spectrum as the default reference.
* Added the option to remove the x axis in an mrs_data plot.
* Added ylim and y_scale options to fit plotting.
* Added %$% operator from magrittr package.
* Added an interpolation option to calc_spec_snr.
* Added hline and vline options to image.mrs_data.

# spant 0.12.0
* Fit results stackplot now has the option to display labels.
* Added the option to reverse eddy current correction.
* Improved GE p-file reader.
* diff function can now be applied to mrs_data objects.
* Complex functions: Re, Im, Mod, Arg and Conj can now be applied to mrs_data 
objects.
* Default simulations for Glc and NAAG have been improved.

# spant 0.11.0
* Added mar argument to plot command.
* td2fd and fd2td now give warnings when used with data already in the target
domain.
* Improved documentation formatting consistency and fixed some spelling errors.
* Added rats method.

# spant 0.10.0
* The names of in-built pulse sequence functions now all start with seq_* to
make them easier to find.
* Added new functions to simulate the following MRS sequences: CPMG, MEGA-PRESS, 
STEAM, sLASER. sLASER sequence kindly contributed by Pierre-Gilles Henry.
* Bug fix for get_mol_names function.
* stackplot function now accepts labels argument and time-domain plotting.
* def_acq_paras function now accepts arguments to override the defaults.
* Added a source field to mol.paras object to cite the origin of the values.
* Option to restore plotting par defaults.
* The magrittr pipe operator is now exported.

# spant 0.9.0
* Updated plotting modes to be one of : "re", "im", "mod" or "arg".
* Updated int_spec function to use "re", "im", or "mod".
* Added a function to replicate data across a particular dimension.
* Added a convenience function to simulate normal looking 1H brain MRS data.
* phase and shift functions now accept vector inputs.

# spant 0.7.0
* Added new function for frequency drift correction.
* Added support for Siemens ima and TWIX SVS data.
* Added support for GE p-file SVS data.
* Added apply_axes fn.
* Support for reading SPM style segmentation results (spm_pve2categorical).

# spant 0.6.0
* Interactive plotting function added for fit results - plot_fit_slice_inter.
* Bug fix for appending dynamic results.
* Bug fix for reading list data files without reference data.
* Bug fix for append_dyns function.
* basis2mrs_data function has been extended to allow the summation of basis
elements and setting of individual amplitudes.
* Added a shift function for manual frequency shift adjustment.
* Added initial unit tests and automated coveralls checking.

# spant 0.5.0
* A default brain PRESS basis is now simulated by the fit_mrs function when the
basis argument isn't specified.
* Added calc_peak_info function for simple singlet analyses.
* crop_spec function now maintains the original frequency scale.
* The basis set used for analyses has now been added to the fit result object.
* Bug fix for simulating basis sets with one element.
* lb function can now be used with basis-set objects.
* Bug fix for spar_sdat reader for non-localised MRS.
* AppVeyor now being used to test Windows compatibility - John Muschelli.

# spant 0.4.0
* Bug fix for SPAR/SDAT SVS voxel dimensions.
* MRSI support added for Philips SPAR/SDAT data.
* Fit plots now default to the full spectral range unless xlim is specified.
* Fit plots allow the x, y, z, coil, dynamic indices to be specified.
* Added the option to subtract the baseline from fit plots.

# spant 0.3.0
* Added stackplot method for fit objects.
* Added functions for registering and visualising SVS volumes on images and 
performing partial volume correction.
* Philips "list data" also now reads noise scans.
* calc_coil_noise_cor, calc_coil_noise_sd functions added to aid coil 
combination.
* Documentation updates for plotting methods.
* Added some simulation methods to userland.

# spant 0.2.0
* Added Siemens RDA format reader.
* Added Philips "list data" format reader.
* Added Bruker paravision format reader.
* Added PROPACK option for HSVD based filtering.
* Added a coil combination function.
* Bug fix for incorrect ppm scale on fit plots when fs != 2000Hz.
* Bug fix for VARPRO analytical jacobian calculation.

# spant 0.1.0
* First public release.
