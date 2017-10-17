# spant 0.6.0
* Interactive plotting function added for fit results - plot_fit_slice_inter.
* Bug fix for appending dynamic results.
* Bug fix for reading list data files without reference data.
* Bug fix for append_dyns function.
* basis2mrs_data function has been extended to allow the summation of basis
elements and setting of individual amplitudes.
* Added a shift function for manual frequency shift adjutment.
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
