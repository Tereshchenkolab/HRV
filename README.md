# HRV
Heart Rate Variability measurement on a single-lead raw Zio Patch ECG signal

## Table of Contents
  - R-peak and Sinus Signal Detection for HRV analysis: R_peak_detector code
  
      Run Main.m file from R_peak_detector folder to calculate r-peaks and select a single 3-minute sinus epoch of signal.
      Utilizes modified pan-tompkins, principle component analysis, and parabolic fitting algorithms for r-peak detection. 
      RR intervals considered sinus if the duration of RR-interval is within 15% of the duration of the previous RR-interval
      
  - HRV Analysis: HRV_full code
  
      Run Main.m file from HRV_full folder on data after running R-peak detection code (above)
      Code will calculate heart rate and HRV parameters (RMSSD, LF, HF, LF/HF, SD1, SD2, SD1/SD2, Sample entropy, renyi entropy)
  
  - Test files: 130812_Z320892913_1563
  
      9 hours of raw data from single lead ZioPatch
      
  - Fully de-identified dataset of the PACE study participants: PACE deidentified ECGpatch HRV dataset 28pts.xls
  
  
### Authors
HRV indices measurement code V.1
Muammar Kabir, PhD, <muammar.kabir@gmail.com>
Nichole Rogovoy, BS <rogovoy@ohsu.edu>
Erick Andres Perez Alday, PhD, <perezald@ohsu.edu>
Annabel Li-Pershing, BS, <lipershi@ohsu.edu>
Larisa Tereshchenko, MD, PhD, <tereshch@ohsu.edu>
  Last update: August 13th, 2019
  
### HRV Indices Calculation
See pdf file in the repository for HRV Indices calculation methods.

### HRV Indices MATLAB Code
See .m files in the repository for HRV parameters calculations. 


### Test files
Test files are provided for HRV calculation testing. Sampling rate 200 Hz. Amplitude resolution 4.88 µV.


### Manuscript to reference:
- Rogovoy NM, Howell SJ, Lee TL, Hamilton C, Perez-Alday EA, Kabir MM, Zhang Y, Kim ED, Fitzpatrick J, Monroy-Trujillo JM, Estrella MM, Sozio SM, Jaar BG, Parekh RS, Tereshchenko LG. Hemodialysis Procedure-Associated Autonomic Imbalance and Cardiac Arrhythmias: Insights From Continuous 14-Day ECG Monitoring. J Am Heart Assoc. 2019 Oct;8(19):e013748. doi: 10.1161/JAHA.119.013748. Epub 2019 Sep 28. PMID: 31564195; PMCID: PMC6806026.  https://www.ncbi.nlm.nih.gov/pubmed/31564195
- BioRxiv preprint: doi: https://doi.org/10.1101/601542
