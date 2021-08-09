---
title: 'spant: An R package for magnetic resonance spectroscopy analysis'
tags:
  - R
  - spectroscopy
  - MRS
  - NMR
  - medical imaging
  - neuroimaging
authors:
  - name: Martin Wilson
    orcid: 0000-0002-2089-3956
    affiliation: 1
affiliations:
 - name: Centre for Human Brain Health and School of Psychology, University of Birmingham, Birmingham, UK
   index: 1
date: 4 August 2021
bibliography: paper.bib
---

# Summary

Magnetic Resonance Spectroscopy (MRS) allows the measurement of small molecules (metabolites) in the body without the use of harmful radiation. Based on the same basic principles and technology behind Magnetic Resonance Imaging (MRI), most modern MRI scanners are also capable of acquiring MRS - making the technique highly suited to a number of clinical applications (@oz:2014). Despite the success of MRS in the research environment, clinical translation has proven slow due to a number of technical and practical reasons, with challenges associated with reliable data processing and analysis having particular importance (@wilson:2019). The `spant` (SPectroscopy ANalysis Tools) package has been developed to: (1) provide open-source implementations of traditional and modern MRS processing and analysis techniques for routine analysis (@near:2021) and (2) aid the development, validation and comparison of new algorithms and analysis pipelines.

# Statement of need

Traditional MRS analysis was dominated by the use of proprietary software, either supplied by scanner manufactures or offline tools such as LCModel (provencher:1993) and jMRUI (naressi:2001). In more recent years there has been a steadily increasing trend toward the use of open-source methods - with some early examples including TARQUIN (@reynolds:2006; @wilson:2011) and AQSES (@poullet:2007). This trend is set to continue with the recent transition of LCModel to an open-source license, and an acceleration in the development of new open-source methods and packages such as FID-A (@simpson:2017), GANNET, VESPA, SUSPECT, OSPRAY and FSL-MRS. The availability of the MRSHub, a new community orientated software sharing and support platform, and the development of the NIfTI MRS file format, to aid data sharing and interoperability, are set to further enhance the ecosystem of open-source MRS analysis tools.

The vast majority of recentely developed open-source MRS analysis tools have been written in either MATLAB or Python. Whilst all languages have strengths and weaknesses, R is particularly suited to the interactive exploration and batch processing of large and complex datasets - typical of MRS and neuroimaging studies. The `spant` package was developed to combine traditional  and modern MRS data processing techniques with strengths of R, including: plotting/visualisation, statistics, machine learning and data wrangling. Furthermore, `spant` may be used to conveniently combine MRS results with other imaging modalities, due to the availabilty of a wide range of R packages focussed on image processing (Neuroconductor ref) and support for the NIfTI data format (clayden, oro.nifti) and its extension to MRS (NIfTI MRS paper).

At the time of writing, `spant` has been used to develop and validate two new MRS spectroscopy analysis algorithms: RATS (ref) and ABfit (ref), and has also been used to study Alzheimer’s Disease (ref) and psychosis (ref) - confirming its suitability for both MRS methods research and clinical studies.

<!---

The availablilty of spant on CRAN provides a straightforward route for installation on Windows, Linux and Mac lowering barr

 - ultimately 

An MRS analysis pipeline is typically composed of the following steps *ref jamie processing paper*:

1) processing : raw signals from the scanner are combined and manipulated to enhance the metabolite signals and supress artefacts. In the case of magnetic resonance spectroscopic imaging (MRSI), spectra are also mapped to spatial locations in a process known as reconstruction.
2) analysis : spectra typically undergo a non-linear fitting procedure to extract unscaled metabolite levels based on prior knowledge of individual metabolite spectra (basis set).
3) quantification : metabolite levels are scaled to meaningful values, such as ratios between two or more metabolites or absolute concentrations.

initatives and tools to support effective data sharing (NIFTI MRS and tools) 

R and RStudio provide a natural environment for this higher level analysis - with advanced and mature tools to interactively organise, visualise and perform a wide range of statistical tests on complex datasets. In addition, reproducible research...

Neuroconductor, ref John Clayden nifti package.

FID/A ref

TARQUIN, LCModel - mainly fitting
OSPRAY, FSL-MRS, VESPA

open methods

emphasis on batch processing

TODO @wilson:2011

Different methods are known to produce different results (Georg paper)

NIfTI MRS format

Data simulation / Monte-Carlo stuff

TODO @simpson:2017

Reproducable reserach - does raw data processing, fitting, MRS visualisation, subject level stat

---!>

# References