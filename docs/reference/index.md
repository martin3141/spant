# Package index

## All functions

- [`Arg(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/Arg.mrs_data.md)
  : Apply Arg operator to an MRS dataset.

- [`Conj(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/Conj.mrs_data.md)
  : Apply Conj operator to an MRS dataset.

- [`Im(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/Im.mrs_data.md)
  : Apply Im operator to an MRS dataset.

- [`Imzap()`](https://martin3141.github.io/spant/reference/Imzap.md) :
  Complex rounding function taken from complexplus package to reduce the
  number of spant dependencies.

- [`Mod(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/Mod.mrs_data.md)
  : Apply Mod operator to an MRS dataset.

- [`Ncoils()`](https://martin3141.github.io/spant/reference/Ncoils.md) :
  Return the total number of coil elements in an MRS dataset.

- [`Ndyns()`](https://martin3141.github.io/spant/reference/Ndyns.md) :
  Return the total number of dynamic scans in an MRS dataset.

- [`Npts()`](https://martin3141.github.io/spant/reference/Npts.md) :
  Return the number of data points in an MRS dataset.

- [`Nspec()`](https://martin3141.github.io/spant/reference/Nspec.md) :
  Return the total number of spectra in an MRS dataset.

- [`Ntrans()`](https://martin3141.github.io/spant/reference/Ntrans.md) :
  Return the total number of acquired transients for an MRS dataset.

- [`Nx()`](https://martin3141.github.io/spant/reference/Nx.md) : Return
  the total number of x locations in an MRS dataset.

- [`Ny()`](https://martin3141.github.io/spant/reference/Ny.md) : Return
  the total number of y locations in an MRS dataset.

- [`Nz()`](https://martin3141.github.io/spant/reference/Nz.md) : Return
  the total number of z locations in an MRS dataset.

- [`Re(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/Re.mrs_data.md)
  : Apply Re operator to an MRS dataset.

- [`abfit_opts()`](https://martin3141.github.io/spant/reference/abfit_opts.md)
  : Return a list of options for an ABfit analysis.

- [`abfit_opts_v1_9_0()`](https://martin3141.github.io/spant/reference/abfit_opts_v1_9_0.md)
  : Return a list of options for an ABfit analysis to maintain
  comparability with analyses performed with version 1.9.0 (and earlier)
  of spant.

- [`abfit_reg_opts()`](https://martin3141.github.io/spant/reference/abfit_reg_opts.md)
  : Return a list of options for an ABfit analysis with regularision.

- [`acquire()`](https://martin3141.github.io/spant/reference/acquire.md)
  : Simulate pulse sequence acquisition.

- [`add_noise()`](https://martin3141.github.io/spant/reference/add_noise.md)
  : Add noise to an mrs_data object.

- [`add_noise_spec_snr()`](https://martin3141.github.io/spant/reference/add_noise_spec_snr.md)
  : Add noise to an mrs_data object to match a given SNR.

- [`align()`](https://martin3141.github.io/spant/reference/align.md) :
  Align spectra to a reference frequency using a convolution based
  method.

- [`apodise_xy()`](https://martin3141.github.io/spant/reference/apodise_xy.md)
  : Apodise MRSI data in the x-y direction with a k-space filter.

- [`append_basis()`](https://martin3141.github.io/spant/reference/append_basis.md)
  : Combine a pair of basis set objects.

- [`append_coils()`](https://martin3141.github.io/spant/reference/append_coils.md)
  : Append MRS data across the coil dimension, assumes they matched
  across the other dimensions.

- [`append_dyns()`](https://martin3141.github.io/spant/reference/append_dyns.md)
  : Append MRS data across the dynamic dimension, assumes they matched
  across the other dimensions.

- [`append_regs()`](https://martin3141.github.io/spant/reference/append_regs.md)
  : Append multiple regressor data frames into a single data frame.

- [`apply_axes()`](https://martin3141.github.io/spant/reference/apply_axes.md)
  : Apply a function over specified array axes.

- [`apply_mrs()`](https://martin3141.github.io/spant/reference/apply_mrs.md)
  : Apply a function across given dimensions of a MRS data object.

- [`apply_pulse()`](https://martin3141.github.io/spant/reference/apply_pulse.md)
  : Simulate an RF pulse on a single spin.

- [`array2mrs_data()`](https://martin3141.github.io/spant/reference/array2mrs_data.md)
  : Convert a 7 dimensional array in into a mrs_data object. The array
  dimensions should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.

- [`assign_dyns()`](https://martin3141.github.io/spant/reference/assign_dyns.md)
  : Overwrite a subset of dynamic scans.

- [`auto_phase()`](https://martin3141.github.io/spant/reference/auto_phase.md)
  : Perform zeroth-order phase correction based on the minimisation of
  the squared difference between the real and magnitude components of
  the spectrum.

- [`auto_phase_bl()`](https://martin3141.github.io/spant/reference/auto_phase_bl.md)
  : Perform zeroth-order phase correction based on expected baseline
  regions.

- [`back_extrap_ar()`](https://martin3141.github.io/spant/reference/back_extrap_ar.md)
  : Back extrapolate time-domain data points using an autoregressive
  model.

- [`basis2dyn_mrs_data()`](https://martin3141.github.io/spant/reference/basis2dyn_mrs_data.md)
  : Convert a basis object to a dynamic mrs_data object.

- [`basis2mrs_data()`](https://martin3141.github.io/spant/reference/basis2mrs_data.md)
  : Convert a basis object to an mrs_data object - where basis signals
  are spread across the dynamic dimension.

- [`bbase()`](https://martin3141.github.io/spant/reference/bbase.md) :
  Generate a spline basis, slightly adapted from : "Splines, knots, and
  penalties", Eilers 2010.

- [`bc_als()`](https://martin3141.github.io/spant/reference/bc_als.md) :
  Baseline correction using the ALS method.

- [`bc_constant()`](https://martin3141.github.io/spant/reference/bc_constant.md)
  : Remove a constant baseline offset based on a reference spectral
  region.

- [`bc_gauss()`](https://martin3141.github.io/spant/reference/bc_gauss.md)
  : Apply and subtract a Gaussian smoother in the spectral domain.

- [`bc_poly()`](https://martin3141.github.io/spant/reference/bc_poly.md)
  : Fit and subtract a polynomial to each spectrum in a dataset.

- [`bc_spline()`](https://martin3141.github.io/spant/reference/bc_spline.md)
  : Fit and subtract a smoothing spline to each spectrum in a dataset.

- [`beta2lw()`](https://martin3141.github.io/spant/reference/beta2lw.md)
  : Covert a beta value in the time-domain to an equivalent linewidth in
  Hz: x \* exp(-i \* t \* t \* beta).

- [`bin_spec()`](https://martin3141.github.io/spant/reference/bin_spec.md)
  : Bin equally spaced spectral regions.

- [`calc_basis_corr_mat()`](https://martin3141.github.io/spant/reference/calc_basis_corr_mat.md)
  : Estimate the correlation matrix for a basis set.

- [`calc_basis_crlbs()`](https://martin3141.github.io/spant/reference/calc_basis_crlbs.md)
  : Estimate the CRLB for each element in a basis set.

- [`calc_coil_noise_cor()`](https://martin3141.github.io/spant/reference/calc_coil_noise_cor.md)
  : Calculate the noise correlation between coil elements.

- [`calc_coil_noise_sd()`](https://martin3141.github.io/spant/reference/calc_coil_noise_sd.md)
  : Calculate the noise standard deviation for each coil element.

- [`calc_design_efficiency()`](https://martin3141.github.io/spant/reference/calc_design_efficiency.md)
  : Calculate the efficiency of a regressor data frame.

- [`calc_ed_from_lambda()`](https://martin3141.github.io/spant/reference/calc_ed_from_lambda.md)
  : Calculate the effective dimensions of a spline smoother from lambda.

- [`calc_peak_info_vec()`](https://martin3141.github.io/spant/reference/calc_peak_info_vec.md)
  : Calculate the FWHM of a peak from a vector of intensity values.

- [`calc_sd_poly()`](https://martin3141.github.io/spant/reference/calc_sd_poly.md)
  : Perform a polynomial fit, subtract and return the standard deviation
  of the residuals.

- [`calc_spec_diff()`](https://martin3141.github.io/spant/reference/calc_spec_diff.md)
  : Calculate the sum of squares differences between two mrs_data
  objects.

- [`calc_spec_snr()`](https://martin3141.github.io/spant/reference/calc_spec_snr.md)
  : Calculate the spectral SNR.

- [`check_lcm()`](https://martin3141.github.io/spant/reference/check_lcm.md)
  : Check LCModel can be run

- [`check_tqn()`](https://martin3141.github.io/spant/reference/check_tqn.md)
  : Check the TARQUIN binary can be run

- [`circ_mask()`](https://martin3141.github.io/spant/reference/circ_mask.md)
  : Create a logical circular mask spanning the full extent of an n x n
  matrix.

- [`coherence_filter()`](https://martin3141.github.io/spant/reference/coherence_filter.md)
  : Zero all coherence orders other than the one supplied as an
  argument.

- [`collapse_to_dyns()`](https://martin3141.github.io/spant/reference/collapse_to_dyns.md)
  : Collapse MRS data by concatenating spectra along the dynamic
  dimension.

- [`comb_coils()`](https://martin3141.github.io/spant/reference/comb_coils.md)
  : Combine coil data based on the first data point of a reference
  signal.

- [`comb_coils_mrsi_gls()`](https://martin3141.github.io/spant/reference/comb_coils_mrsi_gls.md)
  : Combine MRSI coil data using the GLS method presented by An et al
  JMRI 37:1445-1450 (2013).

- [`comb_coils_svs_gls()`](https://martin3141.github.io/spant/reference/comb_coils_svs_gls.md)
  : Combine SVS coil data using the GLS method presented by An et al
  JMRI 37:1445-1450 (2013).

- [`comb_fit_list_dyns()`](https://martin3141.github.io/spant/reference/comb_fit_list_dyns.md)
  : Combine a list of fit results into a single fit result object across
  the dynamic dimension.

- [`comb_fit_list_fit_tables()`](https://martin3141.github.io/spant/reference/comb_fit_list_fit_tables.md)
  : Combine all fitting data points from a list of fits into a single
  data frame.

- [`comb_fit_list_result_tables()`](https://martin3141.github.io/spant/reference/comb_fit_list_result_tables.md)
  : Combine the fit result tables from a list of fit results.

- [`comb_fit_tables()`](https://martin3141.github.io/spant/reference/comb_fit_tables.md)
  : Combine all fitting data points into a single data frame.

- [`comb_metab_ref()`](https://martin3141.github.io/spant/reference/comb_metab_ref.md)
  : Combine a reference and metabolite mrs_data object.

- [`conv_mrs()`](https://martin3141.github.io/spant/reference/conv_mrs.md)
  : Convolve two MRS data objects.

- [`crop_basis()`](https://martin3141.github.io/spant/reference/crop_basis.md)
  :

  Crop `basis_set` object based on a frequency range.

- [`crop_spec()`](https://martin3141.github.io/spant/reference/crop_spec.md)
  :

  Crop `mrs_data` object based on a frequency range.

- [`crop_td_pts()`](https://martin3141.github.io/spant/reference/crop_td_pts.md)
  :

  Crop `mrs_data` object data points in the time-domain.

- [`crop_td_pts_end()`](https://martin3141.github.io/spant/reference/crop_td_pts_end.md)
  :

  Crop `mrs_data` object data points at the end of the FID.

- [`crop_td_pts_pot()`](https://martin3141.github.io/spant/reference/crop_td_pts_pot.md)
  :

  Crop `mrs_data` object data points in the time-domain rounding down to
  the next smallest power of two (pot). Data that already has a pot
  length will not be changed.

- [`crop_xy()`](https://martin3141.github.io/spant/reference/crop_xy.md)
  : Crop an MRSI dataset in the x-y direction

- [`crossprod_3d()`](https://martin3141.github.io/spant/reference/crossprod_3d.md)
  : Compute the vector cross product between vectors x and y. Adapted
  from
  http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function

- [`decimate_mrs_fd()`](https://martin3141.github.io/spant/reference/decimate_mrs_fd.md)
  : Decimate an MRS signal to half the original sampling frequency by
  filtering in the frequency domain before down sampling.

- [`decimate_mrs_td()`](https://martin3141.github.io/spant/reference/decimate_mrs_td.md)
  : Decimate an MRS signal by filtering in the time domain before
  downsampling.

- [`deconv_mrs()`](https://martin3141.github.io/spant/reference/deconv_mrs.md)
  : Deconvolve two MRS data objects.

- [`def_N()`](https://martin3141.github.io/spant/reference/def_N.md) :
  Return the default number of data points in the spectral dimension.

- [`def_acq_paras()`](https://martin3141.github.io/spant/reference/def_acq_paras.md)
  : Return (and optionally modify using the input arguments) a list of
  the default acquisition parameters.

- [`def_fs()`](https://martin3141.github.io/spant/reference/def_fs.md) :
  Return the default sampling frequency in Hz.

- [`def_ft()`](https://martin3141.github.io/spant/reference/def_ft.md) :
  Return the default transmitter frequency in Hz.

- [`def_nuc()`](https://martin3141.github.io/spant/reference/def_nuc.md)
  : Return the default nucleus.

- [`def_ref()`](https://martin3141.github.io/spant/reference/def_ref.md)
  : Return the default reference value for ppm scale.

- [`dicom_reader()`](https://martin3141.github.io/spant/reference/dicom_reader.md)
  : A very simple DICOM reader.

- [`diff_mrs()`](https://martin3141.github.io/spant/reference/diff_mrs.md)
  : Apply the diff operator to an MRS dataset in the FID/spectral
  dimension.

- [`downsample_mrs_fd()`](https://martin3141.github.io/spant/reference/downsample_mrs_fd.md)
  : Downsample an MRS signal by a factor of 2 using an FFT "brick-wall"
  filter.

- [`downsample_mrs_td()`](https://martin3141.github.io/spant/reference/downsample_mrs_td.md)
  : Downsample an MRS signal by a factor of 2 by removing every other
  data point in the time-domain. Note, signals outside the new sampling
  frequency will be aliased.

- [`dyn_acq_times()`](https://martin3141.github.io/spant/reference/dyn_acq_times.md)
  : Return a time scale vector of acquisition times for a dynamic MRS
  scan. The first temporal scan is assigned a value of 0.

- [`ecc()`](https://martin3141.github.io/spant/reference/ecc.md) : Eddy
  current correction.

- [`elliptical_mask()`](https://martin3141.github.io/spant/reference/elliptical_mask.md)
  : Create an elliptical mask stored as a matrix of logical values.

- [`est_noise_sd()`](https://martin3141.github.io/spant/reference/est_noise_sd.md)
  : Estimate the standard deviation of the noise from a segment of an
  mrs_data object.

- [`faceoff()`](https://martin3141.github.io/spant/reference/faceoff.md)
  : Deface a T1 weighted head scan using the FaceOff method described in
  : https://github.com/srikash/FaceOff. Requires ANTs to be installed.

- [`fd2td()`](https://martin3141.github.io/spant/reference/fd2td.md) :
  Transform frequency-domain data to the time-domain.

- [`fd_conv_filt()`](https://martin3141.github.io/spant/reference/fd_conv_filt.md)
  : Frequency-domain convolution based filter.

- [`fd_gauss_smo()`](https://martin3141.github.io/spant/reference/fd_gauss_smo.md)
  : Apply a Gaussian smoother in the spectral domain.

- [`find_bids_mrs()`](https://martin3141.github.io/spant/reference/find_bids_mrs.md)
  : Search for MRS data files in a BIDS filesystem structure.

- [`find_mrs_files()`](https://martin3141.github.io/spant/reference/find_mrs_files.md)
  : Find valid MRS data files recursively from a directory path.

- [`fit_amps()`](https://martin3141.github.io/spant/reference/fit_amps.md)
  :

  Extract the fit amplitudes from an object of class `fit_result`.

- [`fit_asy_pvoigt()`](https://martin3141.github.io/spant/reference/fit_asy_pvoigt.md)
  : Fit a single asymmetric pseudo-Voigt resonance in the frequency
  domain.

- [`fit_diags()`](https://martin3141.github.io/spant/reference/fit_diags.md)
  :

  Calculate diagnostic information for object of class `fit_result`.

- [`fit_mrs()`](https://martin3141.github.io/spant/reference/fit_mrs.md)
  : Perform a fit based analysis of MRS data.

- [`fit_res2csv()`](https://martin3141.github.io/spant/reference/fit_res2csv.md)
  : Write fit results table to a csv file.

- [`fit_svs()`](https://martin3141.github.io/spant/reference/fit_svs.md)
  : Standard SVS 1H brain analysis pipeline.

- [`fit_svs_edited()`](https://martin3141.github.io/spant/reference/fit_svs_edited.md)
  : Edited SVS 1H brain analysis pipeline.

- [`fit_svs_edited_group_results()`](https://martin3141.github.io/spant/reference/fit_svs_edited_group_results.md)
  : Combine edited fitting results for group analysis.

- [`fit_svs_group_results()`](https://martin3141.github.io/spant/reference/fit_svs_group_results.md)
  : Combine fitting results for group analysis.

- [`fit_svs_gui()`](https://martin3141.github.io/spant/reference/fit_svs_gui.md)
  : GUI interface for the standard SVS 1H brain analysis pipeline, this
  is a work in progress, and not ready for serious use.

- [`fit_t1_ti_array()`](https://martin3141.github.io/spant/reference/fit_t1_ti_array.md)
  : Fit a T1 recovery curve, from multiple TIs, to a set of amplitudes.

- [`fit_t1_tr_array()`](https://martin3141.github.io/spant/reference/fit_t1_tr_array.md)
  : Fit a T1 recovery curve, from multiple TRs, to a set of amplitudes.

- [`fit_t2_te_array()`](https://martin3141.github.io/spant/reference/fit_t2_te_array.md)
  : Fit a T2 relaxation curve, from multiple TEs, to a set of
  amplitudes.

- [`fp_phase()`](https://martin3141.github.io/spant/reference/fp_phase.md)
  : Return the phase of the first data point in the time-domain.

- [`fp_phase_correct()`](https://martin3141.github.io/spant/reference/fp_phase_correct.md)
  : Perform a zeroth order phase correction based on the phase of the
  first data point in the time-domain.

- [`fp_scale()`](https://martin3141.github.io/spant/reference/fp_scale.md)
  : Scale the first time-domain data point in an mrs_data object.

- [`fs()`](https://martin3141.github.io/spant/reference/fs.md) : Return
  the sampling frequency in Hz of an MRS dataset.

- [`ft_dyns()`](https://martin3141.github.io/spant/reference/ft_dyns.md)
  : Apply the Fourier transform over the dynamic dimension.

- [`ft_shift()`](https://martin3141.github.io/spant/reference/ft_shift.md)
  : Perform a fft and ffshift on a vector.

- [`ft_shift_mat()`](https://martin3141.github.io/spant/reference/ft_shift_mat.md)
  : Perform a fft and fftshift on a matrix with each column replaced by
  its shifted fft.

- [`gausswin_2d()`](https://martin3141.github.io/spant/reference/gausswin_2d.md)
  : Create a two dimensional Gaussian window function stored as a
  matrix.

- [`gen_F()`](https://martin3141.github.io/spant/reference/gen_F.md) :
  Generate the F product operator.

- [`gen_F_xy()`](https://martin3141.github.io/spant/reference/gen_F_xy.md)
  : Generate the Fxy product operator with a specified phase.

- [`gen_I()`](https://martin3141.github.io/spant/reference/gen_I.md) :
  Generate the I product operator for a single spin.

- [`gen_baseline_reg()`](https://martin3141.github.io/spant/reference/gen_baseline_reg.md)
  : Generate baseline regressor.

- [`gen_bold_reg()`](https://martin3141.github.io/spant/reference/gen_bold_reg.md)
  : Generate BOLD regressors.

- [`gen_conv_reg()`](https://martin3141.github.io/spant/reference/gen_conv_reg.md)
  : Generate regressors by convolving a specified response function with
  a stimulus.

- [`gen_group_reg()`](https://martin3141.github.io/spant/reference/gen_group_reg.md)
  : Expand a regressor matrix for a group analysis.

- [`gen_impulse_reg()`](https://martin3141.github.io/spant/reference/gen_impulse_reg.md)
  : Generate impulse regressors.

- [`gen_numeric_reg()`](https://martin3141.github.io/spant/reference/gen_numeric_reg.md)
  : Generate a regressor from a numeric vector.

- [`gen_poly_reg()`](https://martin3141.github.io/spant/reference/gen_poly_reg.md)
  : Generate polynomial regressors.

- [`gen_trap_reg()`](https://martin3141.github.io/spant/reference/gen_trap_reg.md)
  : Generate trapezoidal regressors.

- [`get_1h_brain_basis_names()`](https://martin3141.github.io/spant/reference/get_1h_brain_basis_names.md)
  : Return a character vector of common 1H molecules found in healthy
  human brain.

- [`get_1h_brain_basis_paras()`](https://martin3141.github.io/spant/reference/get_1h_brain_basis_paras.md)
  :

  Return a list of `mol_parameter` objects suitable for 1H brain MRS
  analyses.

- [`get_1h_brain_basis_paras_v1()`](https://martin3141.github.io/spant/reference/get_1h_brain_basis_paras_v1.md)
  :

  Return a list of `mol_parameter` objects suitable for 1H brain MRS
  analyses.

- [`get_1h_brain_basis_paras_v2()`](https://martin3141.github.io/spant/reference/get_1h_brain_basis_paras_v2.md)
  :

  Return a list of `mol_parameter` objects suitable for 1H brain MRS
  analyses.

- [`get_1h_brain_basis_paras_v3()`](https://martin3141.github.io/spant/reference/get_1h_brain_basis_paras_v3.md)
  :

  Return a list of `mol_parameter` objects suitable for 1H brain MRS
  analyses.

- [`get_1h_braino_basis_names()`](https://martin3141.github.io/spant/reference/get_1h_braino_basis_names.md)
  : Return a character vector of molecules included in the GE BRAINO
  phantom.

- [`get_1h_spectre_basis_names()`](https://martin3141.github.io/spant/reference/get_1h_spectre_basis_names.md)
  : Return a character vector of molecules included in the Gold Star
  Phantoms SPECTRE phantom.

- [`get_2d_psf()`](https://martin3141.github.io/spant/reference/get_2d_psf.md)
  : Get the point spread function (PSF) for a 2D phase encoded MRSI
  scan.

- [`get_acq_paras()`](https://martin3141.github.io/spant/reference/get_acq_paras.md)
  : Return acquisition parameters from a MRS data object.

- [`get_ants_dir()`](https://martin3141.github.io/spant/reference/get_ants_dir.md)
  : Return the ANTs installation directory, or throw an error if not
  found.

- [`get_basis_subset()`](https://martin3141.github.io/spant/reference/get_basis_subset.md)
  : Return a subset of the input basis.

- [`get_dyns()`](https://martin3141.github.io/spant/reference/get_dyns.md)
  : Extract a subset of dynamic scans.

- [`get_even_dyns()`](https://martin3141.github.io/spant/reference/get_even_dyns.md)
  : Return even numbered dynamic scans starting from 1 (2,4,6...).

- [`get_fh_dyns()`](https://martin3141.github.io/spant/reference/get_fh_dyns.md)
  : Return the first half of a dynamic series.

- [`get_fit_map()`](https://martin3141.github.io/spant/reference/get_fit_map.md)
  : Get a data array from a fit result.

- [`get_fp()`](https://martin3141.github.io/spant/reference/get_fp.md) :
  Return the first time-domain data point.

- [`get_guassian_pulse()`](https://martin3141.github.io/spant/reference/get_guassian_pulse.md)
  : Generate a gaussian pulse shape.

- [`get_head_dyns()`](https://martin3141.github.io/spant/reference/get_head_dyns.md)
  : Return the first scans of a dynamic series.

- [`get_hrf()`](https://martin3141.github.io/spant/reference/get_hrf.md)
  : Generate a double gamma model of the HRF as used in SPM.

- [`get_lcm_cmd()`](https://martin3141.github.io/spant/reference/get_lcm_cmd.md)
  : Print the command to run the LCModel command-line program.

- [`get_metab()`](https://martin3141.github.io/spant/reference/get_metab.md)
  : Extract the metabolite component from an mrs_data object.

- [`get_mol_names()`](https://martin3141.github.io/spant/reference/get_mol_names.md)
  :

  Return a character array of names that may be used with the
  `get_mol_paras` function.

- [`get_mol_paras()`](https://martin3141.github.io/spant/reference/get_mol_paras.md)
  :

  Get a `mol_parameters` object for a named molecule.

- [`get_mrs_affine()`](https://martin3141.github.io/spant/reference/get_mrs_affine.md)
  : Generate an affine for nifti generation.

- [`get_mrsi2d_seg()`](https://martin3141.github.io/spant/reference/get_mrsi2d_seg.md)
  : Calculate the partial volume estimates for each voxel in a 2D MRSI
  dataset.

- [`get_mrsi_voi()`](https://martin3141.github.io/spant/reference/get_mrsi_voi.md)
  :

  Generate a MRSI VOI from an `mrs_data` object.

- [`get_mrsi_voxel()`](https://martin3141.github.io/spant/reference/get_mrsi_voxel.md)
  :

  Generate a MRSI voxel from an `mrs_data` object.

- [`get_mrsi_voxel_xy_psf()`](https://martin3141.github.io/spant/reference/get_mrsi_voxel_xy_psf.md)
  :

  Generate a MRSI voxel PSF from an `mrs_data` object.

- [`get_odd_dyns()`](https://martin3141.github.io/spant/reference/get_odd_dyns.md)
  : Return odd numbered dynamic scans starting from 1 (1,3,5...).

- [`get_ref()`](https://martin3141.github.io/spant/reference/get_ref.md)
  : Extract the reference component from an mrs_data object.

- [`get_seg_ind()`](https://martin3141.github.io/spant/reference/get_seg_ind.md)
  : Get the indices of data points lying between two values (end \> x \>
  start).

- [`get_sh_dyns()`](https://martin3141.github.io/spant/reference/get_sh_dyns.md)
  : Return the second half of a dynamic series.

- [`get_slice()`](https://martin3141.github.io/spant/reference/get_slice.md)
  : Return a single slice from a larger MRSI dataset.

- [`get_spant_resources_dir()`](https://martin3141.github.io/spant/reference/get_spant_resources_dir.md)
  : Return the spant resources directory.

- [`get_spin_num()`](https://martin3141.github.io/spant/reference/get_spin_num.md)
  : Return the spin number for a given nucleus.

- [`get_subset()`](https://martin3141.github.io/spant/reference/get_subset.md)
  : Extract a subset of MRS data.

- [`get_svs_voi()`](https://martin3141.github.io/spant/reference/get_svs_voi.md)
  :

  Generate a SVS acquisition volume from an `mrs_data` object.

- [`get_tail_dyns()`](https://martin3141.github.io/spant/reference/get_tail_dyns.md)
  : Return the last scans of a dynamic series.

- [`get_td_amp()`](https://martin3141.github.io/spant/reference/get_td_amp.md)
  : Return an array of amplitudes derived from fitting the initial
  points in the time domain and extrapolating back to t=0.

- [`get_tqn_cmd()`](https://martin3141.github.io/spant/reference/get_tqn_cmd.md)
  : Print the command to run the TARQUIN command-line program.

- [`get_uncoupled_mol()`](https://martin3141.github.io/spant/reference/get_uncoupled_mol.md)
  :

  Generate a `mol_parameters` object for a simple spin system with one
  resonance.

- [`get_voi_cog()`](https://martin3141.github.io/spant/reference/get_voi_cog.md)
  : Calculate the centre of gravity for an image containing 0 and 1's.

- [`get_voi_seg()`](https://martin3141.github.io/spant/reference/get_voi_seg.md)
  : Return the white matter, gray matter and CSF composition of a
  volume.

- [`get_voi_seg_psf()`](https://martin3141.github.io/spant/reference/get_voi_seg_psf.md)
  : Return the white matter, gray matter and CSF composition of a
  volume.

- [`get_voxel()`](https://martin3141.github.io/spant/reference/get_voxel.md)
  : Return a single voxel from a larger mrs dataset.

- [`glm_spec()`](https://martin3141.github.io/spant/reference/glm_spec.md)
  : Perform a GLM analysis of dynamic MRS data in the spectral domain.

- [`glm_spec_fmrs_fl()`](https://martin3141.github.io/spant/reference/glm_spec_fmrs_fl.md)
  : Perform first-level spectral GLM analysis of an fMRS dataset.

- [`glm_spec_fmrs_group()`](https://martin3141.github.io/spant/reference/glm_spec_fmrs_group.md)
  : Perform group-level spectral GLM analysis of an fMRS dataset.

- [`glm_spec_group_linhyp()`](https://martin3141.github.io/spant/reference/glm_spec_group_linhyp.md)
  : Test a group-level spectral GLM linear hypothesis.

- [`grid_shift_xy()`](https://martin3141.github.io/spant/reference/grid_shift_xy.md)
  : Grid shift MRSI data in the x/y dimension.

- [`gridplot()`](https://martin3141.github.io/spant/reference/gridplot.md)
  : Arrange spectral plots in a grid.

- [`gridplot(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/gridplot.mrs_data.md)
  : Arrange spectral plots in a grid.

- [`hsvd()`](https://martin3141.github.io/spant/reference/hsvd.md) :
  HSVD of an mrs_data object.

- [`hsvd_filt()`](https://martin3141.github.io/spant/reference/hsvd_filt.md)
  : HSVD based signal filter.

- [`hsvd_vec()`](https://martin3141.github.io/spant/reference/hsvd_vec.md)
  : HSVD of a complex vector.

- [`hz()`](https://martin3141.github.io/spant/reference/hz.md) : Return
  the frequency scale of an MRS dataset in Hz.

- [`ift_shift()`](https://martin3141.github.io/spant/reference/ift_shift.md)
  : Perform an iffshift and ifft on a vector.

- [`ift_shift_mat()`](https://martin3141.github.io/spant/reference/ift_shift_mat.md)
  : Perform an ifft and ifftshift on a matrix with each column replaced
  by its shifted ifft.

- [`image(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/image.mrs_data.md)
  : Image plot method for objects of class mrs_data.

- [`img2kspace_xy()`](https://martin3141.github.io/spant/reference/img2kspace_xy.md)
  : Transform 2D MRSI data to k-space in the x-y direction.

- [`install_ants()`](https://martin3141.github.io/spant/reference/install_ants.md)
  : Install ANTs / ANTsX from : https://github.com/ANTsX/ANTs/releases

- [`install_cli()`](https://martin3141.github.io/spant/reference/install_cli.md)
  : Install the spant command-line interface scripts to a system path.

- [`install_faceoff()`](https://martin3141.github.io/spant/reference/install_faceoff.md)
  : Install FaceOff from : https://github.com/srikash/FaceOff, requires
  ANTs to be installed to work.

- [`install_oasis_template()`](https://martin3141.github.io/spant/reference/install_oasis_template.md)
  : Install the Oasis brain template from :
  https://doi.org/10.6084/m9.figshare.915436.v2

- [`int_spec()`](https://martin3141.github.io/spant/reference/int_spec.md)
  : Integrate a spectral region.

- [`interleave_dyns()`](https://martin3141.github.io/spant/reference/interleave_dyns.md)
  : Interleave the first and second half of a dynamic series.

- [`inv_even_dyns()`](https://martin3141.github.io/spant/reference/inv_even_dyns.md)
  : Invert even numbered dynamic scans starting from 1 (2,4,6...).

- [`inv_odd_dyns()`](https://martin3141.github.io/spant/reference/inv_odd_dyns.md)
  : Invert odd numbered dynamic scans starting from 1 (1,3,5...).

- [`is.def()`](https://martin3141.github.io/spant/reference/is.def.md) :
  Check if an object is defined, which is the same as being not NULL.

- [`is_fd()`](https://martin3141.github.io/spant/reference/is_fd.md) :
  Check if the chemical shift dimension of an MRS data object is in the
  frequency domain.

- [`kspace2img_xy()`](https://martin3141.github.io/spant/reference/kspace2img_xy.md)
  : Transform 2D MRSI data from k-space to image space in the x-y
  direction.

- [`l2_reg()`](https://martin3141.github.io/spant/reference/l2_reg.md) :
  Perform l2 regularisation artefact suppression.

- [`lb()`](https://martin3141.github.io/spant/reference/lb.md) : Apply
  line-broadening (apodisation) to MRS data or basis object.

- [`lb_renoise()`](https://martin3141.github.io/spant/reference/lb_renoise.md)
  : Apply line-broadening to dynamic MRS data and add normally
  distributed noise (renoise) to reverse any associated improvement in
  SNR.

- [`lofdc()`](https://martin3141.github.io/spant/reference/lofdc.md) :
  Correct linear frequency drift.

- [`lw2alpha()`](https://martin3141.github.io/spant/reference/lw2alpha.md)
  : Covert a linewidth in Hz to an equivalent alpha value in the
  time-domain ie: x \* exp(-t \* alpha).

- [`lw2beta()`](https://martin3141.github.io/spant/reference/lw2beta.md)
  : Covert a linewidth in Hz to an equivalent beta value in the
  time-domain ie: x \* exp(-t \* t \* beta).

- [`make_basis_from_raw()`](https://martin3141.github.io/spant/reference/make_basis_from_raw.md)
  : Make a basis-set object from a directory containing LCModel
  formatted RAW files.

- [`mask_dyns()`](https://martin3141.github.io/spant/reference/mask_dyns.md)
  : Mask an MRS dataset in the dynamic dimension.

- [`mask_fit_res()`](https://martin3141.github.io/spant/reference/mask_fit_res.md)
  : Mask fit result spectra depending on a vector of bool values.

- [`mask_xy()`](https://martin3141.github.io/spant/reference/mask_xy.md)
  : Mask an MRSI dataset in the x-y direction

- [`mask_xy_corners()`](https://martin3141.github.io/spant/reference/mask_xy_corners.md)
  : Mask the four corners of an MRSI dataset in the x-y plane.

- [`mask_xy_ellipse()`](https://martin3141.github.io/spant/reference/mask_xy_ellipse.md)
  : Mask the voxels outside an elliptical region spanning the MRSI
  dataset in the x-y plane.

- [`mask_xy_mat()`](https://martin3141.github.io/spant/reference/mask_xy_mat.md)
  : Mask a 2D MRSI dataset in the x-y dimension.

- [`mat2mrs_data()`](https://martin3141.github.io/spant/reference/mat2mrs_data.md)
  : Convert a matrix (with spectral points in the column dimension and
  dynamics in the row dimensions) into a mrs_data object.

- [`match_files()`](https://martin3141.github.io/spant/reference/match_files.md)
  : Match files based on a vector of input paths and a glob pattern. The
  glob pattern is appended to each path and should match one file only.

- [`match_lineshape()`](https://martin3141.github.io/spant/reference/match_lineshape.md)
  : Apply Voigt line-broadening to match a reference spectrum.

- [`matexp()`](https://martin3141.github.io/spant/reference/matexp.md) :
  Matrix exponential function taken from complexplus package to reduce
  the number of spant dependencies.

- [`max_mrs()`](https://martin3141.github.io/spant/reference/max_mrs.md)
  : Apply the max operator to an MRS dataset.

- [`max_mrs_interp()`](https://martin3141.github.io/spant/reference/max_mrs_interp.md)
  : Apply the max operator to an interpolated MRS dataset.

- [`mean(`*`<list>`*`)`](https://martin3141.github.io/spant/reference/mean.list.md)
  : Calculate the mean spectrum from an mrs_data object.

- [`mean(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/mean.mrs_data.md)
  : Calculate the mean spectrum from an mrs_data object.

- [`mean_dyn_blocks()`](https://martin3141.github.io/spant/reference/mean_dyn_blocks.md)
  : Calculate the mean of adjacent dynamic scans.

- [`mean_dyn_pairs()`](https://martin3141.github.io/spant/reference/mean_dyn_pairs.md)
  : Calculate the pairwise means across a dynamic data set.

- [`mean_dyns()`](https://martin3141.github.io/spant/reference/mean_dyns.md)
  : Calculate the mean dynamic data.

- [`mean_dyns_scheme()`](https://martin3141.github.io/spant/reference/mean_dyns_scheme.md)
  : Average sets of dynamics according to a scheme vector.

- [`mean_mrs_list()`](https://martin3141.github.io/spant/reference/mean_mrs_list.md)
  : Return the mean of a list of mrs_data objects.

- [`mean_vec_blocks()`](https://martin3141.github.io/spant/reference/mean_vec_blocks.md)
  : Calculate the mean of adjacent blocks in a vector.

- [`median_dyns()`](https://martin3141.github.io/spant/reference/median_dyns.md)
  : Calculate the median dynamic data.

- [`mod_td()`](https://martin3141.github.io/spant/reference/mod_td.md) :
  Apply the Modulus operator to the time-domain MRS signal.

- [`mr_data2bids()`](https://martin3141.github.io/spant/reference/mr_data2bids.md)
  : Create a BIDS file structure from a vector of data paths or list of
  mri/mrs data objects.

- [`mrs_data2basis()`](https://martin3141.github.io/spant/reference/mrs_data2basis.md)
  : Convert an mrs_data object to basis object - where basis signals are
  spread across the dynamic dimension in the MRS data.

- [`mrs_data2bids()`](https://martin3141.github.io/spant/reference/mrs_data2bids.md)
  : Create a BIDS file structure from a vector of MRS data paths or list
  of mrs_data objects.

- [`mrs_data2list()`](https://martin3141.github.io/spant/reference/mrs_data2list.md)
  : Split MRS data containing multiple spectra into a list of single
  spectra datasets.

- [`mrs_data2mat()`](https://martin3141.github.io/spant/reference/mrs_data2mat.md)
  : Convert mrs_data object to a matrix, with spectral points in the
  column dimension and dynamics in the row dimension.

- [`mrs_data2spec_mat()`](https://martin3141.github.io/spant/reference/mrs_data2spec_mat.md)
  : Convert mrs_data object to a matrix, with spectral points in the
  column dimension and dynamics in the row dimension.

- [`mrs_data2vec()`](https://martin3141.github.io/spant/reference/mrs_data2vec.md)
  : Convert mrs_data object to a vector.

- [`mvfftshift()`](https://martin3141.github.io/spant/reference/mvfftshift.md)
  : Perform a fftshift on a matrix, with each column replaced by its
  shifted result.

- [`mvifftshift()`](https://martin3141.github.io/spant/reference/mvifftshift.md)
  : Perform an ifftshift on a matrix, with each column replaced by its
  shifted result.

- [`n2coord()`](https://martin3141.github.io/spant/reference/n2coord.md)
  : Print fit coordinates from a single index.

- [`nifti_flip_lr()`](https://martin3141.github.io/spant/reference/nifti_flip_lr.md)
  : Flip the x data dimension order of a nifti image. This corresponds
  to flipping MRI data in the left-right direction, assuming the data in
  save in neurological format (can check with fslorient program).

- [`one_page_pdf()`](https://martin3141.github.io/spant/reference/one_page_pdf.md)
  : Export a one-page pdf of a single fit result

- [`ortho3()`](https://martin3141.github.io/spant/reference/ortho3.md) :
  Display an orthographic projection plot of a nifti object.

- [`ortho3_inter()`](https://martin3141.github.io/spant/reference/ortho3_inter.md)
  : Display an interactive orthographic projection plot of a nifti
  object.

- [`paths2df()`](https://martin3141.github.io/spant/reference/paths2df.md)
  : Split paths into parts based on backslash, forwardslash and dot
  characters and return a data frame of these parts.

- [`peak_info()`](https://martin3141.github.io/spant/reference/peak_info.md)
  : Search for the highest peak in a spectral region and return the
  frequency, height and FWHM.

- [`pg_extrap_xy()`](https://martin3141.github.io/spant/reference/pg_extrap_xy.md)
  : Papoulis-Gerchberg (PG) algorithm method for k-space extrapolation.

- [`phase()`](https://martin3141.github.io/spant/reference/phase.md) :
  Apply phasing parameters to MRS data.

- [`phase_ref_1h_brain()`](https://martin3141.github.io/spant/reference/phase_ref_1h_brain.md)
  : Corrected zero order phase and chemical shift offset in 1H MRS data
  from the brain.

- [`plot(`*`<fit_result>`*`)`](https://martin3141.github.io/spant/reference/plot.fit_result.md)
  :

  Plot the fitting results of an object of class `fit_result`.

- [`plot(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/plot.mrs_data.md)
  : Plotting method for objects of class mrs_data.

- [`plot_bc()`](https://martin3141.github.io/spant/reference/plot_bc.md)
  : Convenience function to plot a baseline estimate with the original
  data.

- [`plot_reg()`](https://martin3141.github.io/spant/reference/plot_reg.md)
  : Plot regressors as an image.

- [`plot_slice_fit()`](https://martin3141.github.io/spant/reference/plot_slice_fit.md)
  : Plot a 2D slice from an MRSI fit result object.

- [`plot_slice_fit_inter()`](https://martin3141.github.io/spant/reference/plot_slice_fit_inter.md)
  : Plot a 2D slice from an MRSI fit result object.

- [`plot_slice_map()`](https://martin3141.github.io/spant/reference/plot_slice_map.md)
  : Plot a slice from a 7 dimensional array.

- [`plot_slice_map_inter()`](https://martin3141.github.io/spant/reference/plot_slice_map_inter.md)
  : Plot an interactive slice map from a data array where voxels can be
  selected to display a corresponding spectrum.

- [`plot_spec_sd()`](https://martin3141.github.io/spant/reference/plot_spec_sd.md)
  : Plot the spectral standard deviation.

- [`plot_voi_overlay()`](https://martin3141.github.io/spant/reference/plot_voi_overlay.md)
  : Plot a volume as an image overlay.

- [`plot_voi_overlay_seg()`](https://martin3141.github.io/spant/reference/plot_voi_overlay_seg.md)
  : Plot a volume as an overlay on a segmented brain volume.

- [`ppm()`](https://martin3141.github.io/spant/reference/ppm.md) :
  Return the ppm scale of an MRS dataset or fit result.

- [`precomp()`](https://martin3141.github.io/spant/reference/precomp.md)
  : Save function results to file and load on subsequent calls to avoid
  repeat computation.

- [`preproc_svs()`](https://martin3141.github.io/spant/reference/preproc_svs.md)
  : Preprocess and perform quality assessment of a single SVS data set.

- [`preproc_svs_dataset()`](https://martin3141.github.io/spant/reference/preproc_svs_dataset.md)
  : Preprocess and perform quality assessment of one or more SVS data
  sets.

- [`print(`*`<fit_result>`*`)`](https://martin3141.github.io/spant/reference/print.fit_result.md)
  :

  Print a summary of an object of class `fit_result`.

- [`print(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/print.mrs_data.md)
  : Print a summary of mrs_data parameters.

- [`qn_states()`](https://martin3141.github.io/spant/reference/qn_states.md)
  : Get the quantum coherence matrix for a spin system.

- [`rats()`](https://martin3141.github.io/spant/reference/rats.md) :
  Robust Alignment to a Target Spectrum (RATS).

- [`re_weighting()`](https://martin3141.github.io/spant/reference/re_weighting.md)
  : Apply a weighting to the FID to enhance spectral resolution.

- [`read_basis()`](https://martin3141.github.io/spant/reference/read_basis.md)
  : Read a basis file in LCModel .basis format.

- [`read_basis_niidir()`](https://martin3141.github.io/spant/reference/read_basis_niidir.md)
  : Read a basis folder containing one NIfTI MRS file for each basis
  element.

- [`read_dkd_moco_log()`](https://martin3141.github.io/spant/reference/read_dkd_moco_log.md)
  : Read the log from Dinesh's MoCo sLASER sequence.

- [`read_ima_coil_dir()`](https://martin3141.github.io/spant/reference/read_ima_coil_dir.md)
  : Read a directory containing Siemens MRS IMA files and combine along
  the coil dimension. Note that the coil ID is inferred from the sorted
  file name and should be checked when consistency is required between
  two directories.

- [`read_ima_dyn_dir()`](https://martin3141.github.io/spant/reference/read_ima_dyn_dir.md)
  : Read a directory containing Siemens MRS IMA files and combine along
  the dynamic dimension. Note that the coil ID is inferred from the
  sorted file name and should be checked when consistency is required.

- [`read_lcm_coord()`](https://martin3141.github.io/spant/reference/read_lcm_coord.md)
  : Read an LCModel formatted coord file containing fit information.

- [`read_mrs()`](https://martin3141.github.io/spant/reference/read_mrs.md)
  : Read MRS data from the filesystem.

- [`read_mrs_tqn()`](https://martin3141.github.io/spant/reference/read_mrs_tqn.md)
  : Read MRS data using the TARQUIN software package.

- [`read_pulse_ascii()`](https://martin3141.github.io/spant/reference/read_pulse_ascii.md)
  : Read an ASCII formatted pulse file.

- [`read_pulse_bruker()`](https://martin3141.github.io/spant/reference/read_pulse_bruker.md)
  : Read a Bruker formatted pulse file

- [`read_pulse_pta()`](https://martin3141.github.io/spant/reference/read_pulse_pta.md)
  : Read a .pta formatted pulse file compatible with Siemens PulseTool.

- [`read_siemens_txt_hdr()`](https://martin3141.github.io/spant/reference/read_siemens_txt_hdr.md)
  : Read the text format header found in Siemens IMA and TWIX data
  files.

- [`read_tqn_fit()`](https://martin3141.github.io/spant/reference/read_tqn_fit.md)
  : Reader for csv fit results generated by TARQUIN.

- [`read_tqn_result()`](https://martin3141.github.io/spant/reference/read_tqn_result.md)
  : Reader for csv results generated by TARQUIN.

- [`recon_imag()`](https://martin3141.github.io/spant/reference/recon_imag.md)
  : Reconstruct complex time-domain data from the real part of
  frequency-domain data.

- [`recon_imag_vec()`](https://martin3141.github.io/spant/reference/recon_imag_vec.md)
  : Reconstruct complex time-domain data from the real part of
  frequency-domain data.

- [`recon_twix_2d_mrsi()`](https://martin3141.github.io/spant/reference/recon_twix_2d_mrsi.md)
  : Reconstruct 2D MRSI data from a twix file loaded with read_mrs.

- [`rectangular_mask()`](https://martin3141.github.io/spant/reference/rectangular_mask.md)
  : Create a rectangular mask stored as a matrix of logical values.

- [`rep_array_dim()`](https://martin3141.github.io/spant/reference/rep_array_dim.md)
  : Repeat an array over a given dimension.

- [`rep_dyn()`](https://martin3141.github.io/spant/reference/rep_dyn.md)
  : Replicate a scan in the dynamic dimension.

- [`rep_mrs()`](https://martin3141.github.io/spant/reference/rep_mrs.md)
  : Replicate a scan over a given dimension.

- [`resample_basis()`](https://martin3141.github.io/spant/reference/resample_basis.md)
  : Resample a basis-set to match a mrs_data acquisition.

- [`resample_img()`](https://martin3141.github.io/spant/reference/resample_img.md)
  : Resample an image to match a target image space.

- [`resample_voi()`](https://martin3141.github.io/spant/reference/resample_voi.md)
  : Resample a VOI to match a target image space using nearest-neighbour
  interpolation.

- [`reslice_to_mrs()`](https://martin3141.github.io/spant/reference/reslice_to_mrs.md)
  : Reslice a nifti object to match the orientation of mrs data.

- [`reson_table2mrs_data()`](https://martin3141.github.io/spant/reference/reson_table2mrs_data.md)
  : Generate mrs_data from a table of single Lorentzian resonances.

- [`rm_dyns()`](https://martin3141.github.io/spant/reference/rm_dyns.md)
  : Remove a subset of dynamic scans.

- [`scale_amp_legacy()`](https://martin3141.github.io/spant/reference/scale_amp_legacy.md)
  : Apply water reference scaling to a fitting results object to yield
  metabolite quantities in units of "mmol per Kg wet weight".

- [`scale_amp_molal()`](https://martin3141.github.io/spant/reference/scale_amp_molal.md)
  : Apply water reference scaling to a fitting results object to yield
  metabolite quantities in millimolar (mM) units (mol / kg of tissue
  water).

- [`scale_amp_molal_pvc()`](https://martin3141.github.io/spant/reference/scale_amp_molal_pvc.md)
  : Apply water reference scaling to a fitting results object to yield
  metabolite quantities in millimolar (mM) units (mol / kg of tissue
  water).

- [`scale_amp_molar()`](https://martin3141.github.io/spant/reference/scale_amp_molar.md)
  : Apply water reference scaling to a fitting results object to yield
  metabolite quantities in millimolar (mM) units (mol / Litre of
  tissue). This function is depreciated, please use scale_amp_legacy
  instead.

- [`scale_amp_molar2molal_pvc()`](https://martin3141.github.io/spant/reference/scale_amp_molar2molal_pvc.md)
  : Convert default LCM/TARQUIN concentration scaling to molal units
  with partial volume correction.

- [`scale_amp_ratio()`](https://martin3141.github.io/spant/reference/scale_amp_ratio.md)
  : Scale fitted amplitudes to a ratio of signal amplitude.

- [`scale_amp_ratio_value()`](https://martin3141.github.io/spant/reference/scale_amp_ratio_value.md)
  : Scale fitted amplitudes to a ratio of signal amplitude.

- [`scale_amp_water_ratio()`](https://martin3141.github.io/spant/reference/scale_amp_water_ratio.md)
  : Scale metabolite amplitudes as a ratio to the unsuppressed water
  amplitude.

- [`scale_basis_amp()`](https://martin3141.github.io/spant/reference/scale_basis_amp.md)
  : Scale a basis object by a scalar.

- [`scale_basis_from_singlet()`](https://martin3141.github.io/spant/reference/scale_basis_from_singlet.md)
  : Scale a basis-set to be consistent with spant assumptions for water
  scaling.

- [`scale_mrs_amp()`](https://martin3141.github.io/spant/reference/scale_mrs_amp.md)
  : Scale an mrs_data object by a scalar or vector or amplitudes.

- [`scale_spec()`](https://martin3141.github.io/spant/reference/scale_spec.md)
  : Scale mrs_data to a spectral region.

- [`sd()`](https://martin3141.github.io/spant/reference/sd.md) :
  Calculate the standard deviation spectrum from an mrs_data object.

- [`sd(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/sd.mrs_data.md)
  : Calculate the standard deviation spectrum from an mrs_data object.

- [`seconds()`](https://martin3141.github.io/spant/reference/seconds.md)
  : Return a time scale vector to match the FID of an MRS data object.

- [`segment_t1_ants()`](https://martin3141.github.io/spant/reference/segment_t1_ants.md)
  : Segment T1 weighted MRI data using the ANTs binaries and write to
  file.

- [`segment_t1_fsl()`](https://martin3141.github.io/spant/reference/segment_t1_fsl.md)
  : Segment T1 weighted MRI data using FSL FAST and write to file. Runs
  bet as a preprocessing step by default.

- [`segment_t1_rpyants()`](https://martin3141.github.io/spant/reference/segment_t1_rpyants.md)
  : Segment T1 weighted MRI data using the rpyANTs interface to ANTs and
  write to file.

- [`seq_cpmg_ideal()`](https://martin3141.github.io/spant/reference/seq_cpmg_ideal.md)
  : CPMG style sequence with ideal pulses.

- [`seq_mega_press_ideal()`](https://martin3141.github.io/spant/reference/seq_mega_press_ideal.md)
  : MEGA-PRESS sequence with ideal localisation pulses and Gaussian
  shaped editing pulse.

- [`seq_press_2d_shaped()`](https://martin3141.github.io/spant/reference/seq_press_2d_shaped.md)
  : PRESS sequence with shaped refocusing pulses.

- [`seq_press_ideal()`](https://martin3141.github.io/spant/reference/seq_press_ideal.md)
  : PRESS sequence with ideal pulses.

- [`seq_pulse_acquire()`](https://martin3141.github.io/spant/reference/seq_pulse_acquire.md)
  : Simple pulse and acquire sequence with ideal pulses.

- [`seq_slaser_ideal()`](https://martin3141.github.io/spant/reference/seq_slaser_ideal.md)
  : sLASER sequence with ideal pulses.

- [`seq_spin_echo_ideal()`](https://martin3141.github.io/spant/reference/seq_spin_echo_ideal.md)
  : Spin echo sequence with ideal pulses.

- [`seq_steam_ideal()`](https://martin3141.github.io/spant/reference/seq_steam_ideal.md)
  : STEAM sequence with ideal pulses.

- [`seq_steam_ideal_cof()`](https://martin3141.github.io/spant/reference/seq_steam_ideal_cof.md)
  : STEAM sequence with ideal pulses and coherence order filtering to
  simulate gradient crushers.

- [`seq_steam_ideal_young()`](https://martin3141.github.io/spant/reference/seq_steam_ideal_young.md)
  : STEAM sequence with ideal pulses using the z-rotation gradient
  simulation method described by Young et al JMR 140, 146-152 (1999).

- [`set_Ntrans()`](https://martin3141.github.io/spant/reference/set_Ntrans.md)
  : Set the number of transients for an mrs_data object.

- [`set_ants_dir()`](https://martin3141.github.io/spant/reference/set_ants_dir.md)
  : Set the ANTs installation directory location.

- [`set_def_acq_paras()`](https://martin3141.github.io/spant/reference/set_def_acq_paras.md)
  : Set the default acquisition parameters.

- [`set_lcm_cmd()`](https://martin3141.github.io/spant/reference/set_lcm_cmd.md)
  : Set the command to run the LCModel command-line program.

- [`set_lw()`](https://martin3141.github.io/spant/reference/set_lw.md) :
  Apply line-broadening to an mrs_data object to achieve a specified
  linewidth.

- [`set_mask_xy_mat()`](https://martin3141.github.io/spant/reference/set_mask_xy_mat.md)
  : Set the masked voxels in a 2D MRSI dataset to given spectrum.

- [`set_precomp_mode()`](https://martin3141.github.io/spant/reference/set_precomp_mode.md)
  : Set the precompute mode.

- [`set_precomp_verbose()`](https://martin3141.github.io/spant/reference/set_precomp_verbose.md)
  : Set the verbosity of the precompute function.

- [`set_ref()`](https://martin3141.github.io/spant/reference/set_ref.md)
  : Set the ppm reference value (eg ppm value at 0Hz).

- [`set_td_pts()`](https://martin3141.github.io/spant/reference/set_td_pts.md)
  : Set the number of time-domain data points, truncating or
  zero-filling as appropriate.

- [`set_tqn_cmd()`](https://martin3141.github.io/spant/reference/set_tqn_cmd.md)
  : Set the command to run the TARQUIN command-line program.

- [`set_tr()`](https://martin3141.github.io/spant/reference/set_tr.md) :
  Set the repetition time of an MRS dataset.

- [`shift()`](https://martin3141.github.io/spant/reference/shift.md) :
  Apply a frequency shift to MRS data.

- [`shift_basis()`](https://martin3141.github.io/spant/reference/shift_basis.md)
  : Apply frequency shifts to basis set signals.

- [`sim_asy_pvoigt()`](https://martin3141.github.io/spant/reference/sim_asy_pvoigt.md)
  : Generate an asymmetric pseudo-Voigt resonance in the frequency
  domain.

- [`sim_basis()`](https://martin3141.github.io/spant/reference/sim_basis.md)
  : Simulate a basis set object.

- [`sim_basis_1h_brain()`](https://martin3141.github.io/spant/reference/sim_basis_1h_brain.md)
  : Simulate a basis-set suitable for 1H brain MRS analysis acquired
  with a PRESS sequence. Note, ideal pulses are assumed.

- [`sim_basis_1h_brain_press()`](https://martin3141.github.io/spant/reference/sim_basis_1h_brain_press.md)
  : Simulate a basis-set suitable for 1H brain MRS analysis acquired
  with a PRESS sequence. Note, ideal pulses are assumed.

- [`sim_basis_mm_lip_lcm()`](https://martin3141.github.io/spant/reference/sim_basis_mm_lip_lcm.md)
  : Simulate a macromolecular and lipid basis-set suitable for 1H brain
  MRS analysis.

- [`sim_basis_tqn()`](https://martin3141.github.io/spant/reference/sim_basis_tqn.md)
  : Simulate a basis file using TARQUIN.

- [`sim_brain_1h()`](https://martin3141.github.io/spant/reference/sim_brain_1h.md)
  : Simulate MRS data with a similar appearance to normal brain (by
  default).

- [`sim_mol()`](https://martin3141.github.io/spant/reference/sim_mol.md)
  :

  Simulate a `mol_parameter` object.

- [`sim_noise()`](https://martin3141.github.io/spant/reference/sim_noise.md)
  : Simulate an mrs_data object containing simulated Gaussian noise.

- [`sim_resonances()`](https://martin3141.github.io/spant/reference/sim_resonances.md)
  : Simulate a MRS data object containing a set of simulated resonances.

- [`sim_th_excit_profile()`](https://martin3141.github.io/spant/reference/sim_th_excit_profile.md)
  : Simulate an ideal pulse excitation profile by smoothing a top-hat
  function with a Gaussian.

- [`sim_zero()`](https://martin3141.github.io/spant/reference/sim_zero.md)
  : Simulate an mrs_data object containing complex zero valued samples.

- [`smooth_dyns()`](https://martin3141.github.io/spant/reference/smooth_dyns.md)
  : Smooth data across the dynamic dimension with a Gaussian kernel.

- [`sort_basis()`](https://martin3141.github.io/spant/reference/sort_basis.md)
  : Sort the basis-set elements alphabetically.

- [`spant`](https://martin3141.github.io/spant/reference/spant-package.md)
  [`spant-package`](https://martin3141.github.io/spant/reference/spant-package.md)
  : spant: spectroscopy analysis tools.

- [`spant_abfit_benchmark()`](https://martin3141.github.io/spant/reference/spant_abfit_benchmark.md)
  : Simulate and fit some spectra with ABfit for benchmarking purposes.
  Basic timing and performance metrics will be printed.

- [`spant_sim_fmrs_dataset()`](https://martin3141.github.io/spant/reference/spant_sim_fmrs_dataset.md)
  : Simulate an example fMRS dataset for a block design fMRS experiment
  and export a BIDS structure.

- [`spant_simulation_benchmark()`](https://martin3141.github.io/spant/reference/spant_simulation_benchmark.md)
  : Simulate a typical metabolite basis set for benchmarking. Timing
  metrics will be printed on completion.

- [`spec_decomp()`](https://martin3141.github.io/spant/reference/spec_decomp.md)
  : Decompose an mrs_data object into white and gray matter spectra.

- [`spec_op()`](https://martin3141.github.io/spant/reference/spec_op.md)
  : Perform a mathematical operation on a spectral region.

- [`spin_sys()`](https://martin3141.github.io/spant/reference/spin_sys.md)
  : Create a spin system object for pulse sequence simulation.

- [`spm_pve2categorical()`](https://martin3141.github.io/spant/reference/spm_pve2categorical.md)
  : Convert SPM style segmentation files to a single categorical image
  where the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.

- [`ssp()`](https://martin3141.github.io/spant/reference/ssp.md) :
  Signal space projection method for lipid suppression.

- [`stackplot()`](https://martin3141.github.io/spant/reference/stackplot.md)
  : Produce a plot with multiple traces.

- [`stackplot(`*`<fit_result>`*`)`](https://martin3141.github.io/spant/reference/stackplot.fit_result.md)
  :

  Plot the fitting results of an object of class `fit_result` with
  individual basis set components shown.

- [`stackplot(`*`<mrs_data>`*`)`](https://martin3141.github.io/spant/reference/stackplot.mrs_data.md)
  : Stackplot plotting method for objects of class mrs_data.

- [`sub_even_from_odd_dyns()`](https://martin3141.github.io/spant/reference/sub_even_from_odd_dyns.md)
  : Subtract the even dynamic scans from the odd dynamic scans.

- [`sub_first_dyn()`](https://martin3141.github.io/spant/reference/sub_first_dyn.md)
  : Subtract the first dynamic spectrum from a dynamic series.

- [`sub_mean_dyns()`](https://martin3141.github.io/spant/reference/sub_mean_dyns.md)
  : Subtract the mean dynamic spectrum from a dynamic series.

- [`sub_median_dyns()`](https://martin3141.github.io/spant/reference/sub_median_dyns.md)
  : Subtract the median dynamic spectrum from a dynamic series.

- [`sub_odd_from_even_dyns()`](https://martin3141.github.io/spant/reference/sub_odd_from_even_dyns.md)
  : Subtract the odd dynamic scans from the even dynamic scans.

- [`subtract_rest_task()`](https://martin3141.github.io/spant/reference/subtract_rest_task.md)
  : Subtract mean rest spectrum from mean task spectrum after applying
  optimal linebroadening to the mean task spectrum. Usually used to
  correct for the BOLD lineshape narrowing effect in 1H fMRS data.

- [`sum_coils()`](https://martin3141.github.io/spant/reference/sum_coils.md)
  : Calculate the sum across receiver coil elements.

- [`sum_dyns()`](https://martin3141.github.io/spant/reference/sum_dyns.md)
  : Calculate the sum of data dynamics.

- [`sum_mrs()`](https://martin3141.github.io/spant/reference/sum_mrs.md)
  : Sum two mrs_data objects.

- [`sum_mrs_list()`](https://martin3141.github.io/spant/reference/sum_mrs_list.md)
  : Return the sum of a list of mrs_data objects.

- [`sv_res_table()`](https://martin3141.github.io/spant/reference/sv_res_table.md)
  : Output a table of fit amplitudes and error estimates for a
  single-voxel fit.

- [`svs_1h_brain_analysis()`](https://martin3141.github.io/spant/reference/svs_1h_brain_analysis.md)
  : Standard SVS 1H brain analysis pipeline.

- [`svs_1h_brain_analysis_dev()`](https://martin3141.github.io/spant/reference/svs_1h_brain_analysis_dev.md)
  : Standard SVS 1H brain analysis pipeline.

- [`svs_1h_brain_batch_analysis()`](https://martin3141.github.io/spant/reference/svs_1h_brain_batch_analysis.md)
  : Batch interface to the standard SVS 1H brain analysis pipeline.

- [`t_test_spec()`](https://martin3141.github.io/spant/reference/t_test_spec.md)
  : Perform a t-test on spectral data points.

- [`td2fd()`](https://martin3141.github.io/spant/reference/td2fd.md) :
  Transform time-domain data to the frequency-domain.

- [`td_conv_filt()`](https://martin3141.github.io/spant/reference/td_conv_filt.md)
  : Time-domain convolution based filter.

- [`tdsr()`](https://martin3141.github.io/spant/reference/tdsr.md) :
  Time-domain spectral registration.

- [`te()`](https://martin3141.github.io/spant/reference/te.md) : Return
  the echo time of an MRS dataset.

- [`tr()`](https://martin3141.github.io/spant/reference/tr.md) : Return
  the repetition time of an MRS dataset.

- [`trim_paths()`](https://martin3141.github.io/spant/reference/trim_paths.md)
  : Trim a vector of filesystem paths.

- [`varpro_3_para_opts()`](https://martin3141.github.io/spant/reference/varpro_3_para_opts.md)
  : Return a list of options for VARPRO based fitting with 3 free
  parameters.

- [`varpro_basic_opts()`](https://martin3141.github.io/spant/reference/varpro_basic_opts.md)
  : Return a list of options for a basic VARPRO analysis.

- [`varpro_opts()`](https://martin3141.github.io/spant/reference/varpro_opts.md)
  : Return a list of options for VARPRO based fitting.

- [`vec2mrs_data()`](https://martin3141.github.io/spant/reference/vec2mrs_data.md)
  : Convert a vector into a mrs_data object.

- [`write_basis()`](https://martin3141.github.io/spant/reference/write_basis.md)
  : Write a basis object to an LCModel .basis formatted file.

- [`write_basis_niidir()`](https://martin3141.github.io/spant/reference/write_basis_niidir.md)
  : Write a basis object to folder containing one NIfTI MRS file for
  each basis element.

- [`write_basis_tqn()`](https://martin3141.github.io/spant/reference/write_basis_tqn.md)
  : Generate a basis file using TARQUIN.

- [`write_mrs()`](https://martin3141.github.io/spant/reference/write_mrs.md)
  : Write MRS data object to file.

- [`write_mrs_nifti()`](https://martin3141.github.io/spant/reference/write_mrs_nifti.md)
  : Write MRS data object to file in NIFTI format.

- [`write_pulse_ascii()`](https://martin3141.github.io/spant/reference/write_pulse_ascii.md)
  : Write an ASCII formatted pulse file.

- [`zero_fade_spec()`](https://martin3141.github.io/spant/reference/zero_fade_spec.md)
  : Fade a spectrum to zero by frequency domain multiplication with a
  tanh function. Note this operation distorts data points at the end of
  the FID.

- [`zero_higher_orders()`](https://martin3141.github.io/spant/reference/zero_higher_orders.md)
  : Zero all coherences including and above a given order.

- [`zero_td_pts_end()`](https://martin3141.github.io/spant/reference/zero_td_pts_end.md)
  :

  Set `mrs_data` object data points at the end of the FID to zero.

- [`zf()`](https://martin3141.github.io/spant/reference/zf.md) :
  Zero-fill MRS data in the time domain.

- [`zf_xy()`](https://martin3141.github.io/spant/reference/zf_xy.md) :
  Zero-fill MRSI data in the k-space x-y direction.
