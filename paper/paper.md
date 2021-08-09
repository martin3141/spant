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

Magnetic Resonance Spectroscopy (MRS) allows the measurement of small molecules (metabolites) in the body without the use of harmful radiation. Based on the same basic principles and technology behind Magnetic Resonance Imaging (MRI), most modern MRI scanners are also capable of acquiring MRS — making the technique highly suited to a number of clinical applications (@oz:2014). Despite the success of MRS in the research environment, clinical translation has proven slow due to a number of technical and practical reasons, with challenges associated with reliable data processing and analysis having particular importance (@wilson:2019a). The `spant` (SPectroscopy ANalysis Tools) package has been developed to: (1) provide open-source implementations of traditional and modern MRS processing and analysis techniques for routine analysis (@near:2021) and (2) aid the development, validation and comparison of new algorithms and analysis pipelines.

# Statement of need

Traditional MRS analysis was dominated by the use of proprietary software, either supplied by scanner manufactures or offline tools such as LCModel (@provencher:1993) and jMRUI (@naressi:2001). In more recent years there has been a steadily increasing trend toward the use of open-source methods — with some early examples including TARQUIN (@reynolds:2006; @wilson:2011) and AQSES (@poullet:2007). This trend is set to continue with the recent transition of LCModel to an open-source license, and an acceleration in the development of new open-source methods and packages such as Vespa (@soher:2011), Gannet (@edden:2014), FID-A (@simpson:2017), Osprey (@oeltzschner:2020), suspect (@rowland:2021) and FSL-MRS (@clarke:2021a). The availability of the MRSHub <https://mrshub.org/>, a new community orientated software sharing and support platform, and the recent development of the NIfTI MRS file format (@clarke:2021b), to aid data sharing and interoperability, are set to further enhance the ecosystem of open-source MRS analysis tools.

The vast majority of recentely developed open-source MRS analysis tools have been written in either MATLAB or Python. Whilst all languages have strengths and weaknesses, R is particularly suited to the interactive exploration and batch processing of large and complex datasets — typical of MRS and neuroimaging studies. The `spant` package was developed to combine traditional  and modern MRS data processing techniques with strengths of R, including: plotting/visualisation, statistics, machine learning and data wrangling. Furthermore, `spant` may be used to conveniently combine MRS results with other imaging modalities, due to the availabilty of a wide range of R packages focussed on image processing (@muschelli:2019) and support for the NIfTI data format (@whitcher:2011; @clayden:2021).

At the time of writing, `spant` has been used to develop and validate two new MRS spectroscopy analysis algorithms: RATS (@wilson:2019b) and ABfit (@wilson:2021), and has also been used to study cancer (@franco:2021), Alzheimer’s Disease (@montal:2021) and psychosis (@fisher:2020) — confirming its suitability for both MRS methods research and clinical studies.

# Acknowledgements

Particular thanks go to Dr Jonathan D. Clayden and Dr Robert W. Cox for their work on the `RNifti` package (@clayden:2021) and NIfTI standard <https://nifti.nimh.nih.gov/> — which have substantially expanded the capabilites of `spant`.

# References