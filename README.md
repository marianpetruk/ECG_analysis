# ECG ANALYSIS

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)


## Motivation
ECG is a type of time-series data, it has its own particular properties.
This project is created to learn these specifics and how to obtain valuable features from the ECG signal. In particular, from QRS complexes, RR intervals.


## QRS complex
Visualization:

<img src="images/ECG_principle_slow.gif" width="300px"/>


## Screenshots

Raw signal:
![Raw signal](images/raw_signal.png)

Filter frequency response:
![Filter frequency response](images/Filter_frequency_response.png)

Applied filter:
![Applied filter](images/applied_filter.png)

R peaks detection process:
![R peaks detection process](images/R_peaks_detection_process.png)

Comparison with different R-peaks detecors:
![Comparison with different R-peaks detecors](images/Comparison_with_different_R-peaks_detecors.png)




## Dependencies
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/nb_conda/badges/installer/conda.svg)](https://anaconda.org/anaconda/nb_conda): `conda install nb_conda`
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/matplotlib/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda): `conda install matplotlib`
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/seaborn/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda): `conda install seaborn`
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/numpy/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda): `conda install numpy`
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/pywavelets/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda): `conda install pywavelets`
  - [![Anaconda-Server Badge](https://anaconda.org/anaconda/scipy/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda): `conda install scipy`
  - [![PyPI version](https://badge.fury.io/py/biosppy.svg)](https://badge.fury.io/py/biosppy): `pip install biosppy`
  - [![PyPI version](https://badge.fury.io/py/pyhrv.svg)](https://badge.fury.io/py/pyhrv): `pip install pyhrv`

## Literature:
1. [![DOI:10.1109/TBME.1985.325532](https://zenodo.org/badge/DOI/10.1109/TBME.1985.325532.svg)](https://doi.org/10.1109/TBME.1985.325532) Pan-Tomkins algorithm (Pan J., Tompkins W. J., A real-time QRS detection algorithm, IEEE Transactions on Biomedical Engineering, Vol. BME-32, No. 3, March 1985, pp. 230-236).
2. [ECGwaves.com](https://ecgwaves.com/ecg-normal-p-wave-qrs-complex-st-segment-t-wave-j-point/)
3. "ECG filtering T-61.181" – Biomedical Signal Processing Presentation 11.11.2004 Matti Aksela (Aalto University)
4.  "Mining the ECG: Algorithms and Applications" - 2015 KU Leuven – Faculty of Engineering Science Carolina Varon
5. https://imotions.com/blog/heart-rate-variability/
6. http://www.medteq.info/med/ECGFilters
7. http://www.ems12lead.com/2014/03/10/understanding-ecg-filtering/
8. [Band-pass filter | Wikipedia](https://en.wikipedia.org/wiki/Butterworth_filter)
9. [Butterworth filter | Wikipedia](https://en.wikipedia.org/wiki/Butterworth_filter)
10. [Heart rate variability (HRV) | Wikipedia](https://en.wikipedia.org/wiki/Heart_rate_variability)
11. [Electrocardiography | Wikipedia](https://en.wikipedia.org/wiki/Electrocardiography)
12.  [![DOI:10.1093/acprof:oso/9780195058239.003.0019](https://zenodo.org/badge/DOI/10.1093/acprof:oso/9780195058239.003.0019.svg)](https://doi.org/10.1093/acprof:oso/9780195058239.003.0019) Bioelectromagnetism. 19. The Basis of ECG Diagnosis 1995 by Jaakko Malmivuo, Robert Plonsey
13. [Exploring Heart Rate Variability using Python](https://blog.orikami.nl/exploring-heart-rate-variability-using-python-483a7037c64d)
14. [![DOI:10.1037/1089-2680.10.3.229](https://zenodo.org/badge/DOI/10.1037/1089-2680.10.3.229.svg)](https://doi.org/10.1037/1089-2680.10.3.229) Heart rate variability as an index of regulated emotional responding (Appelhans, Bradley M.,Luecken, Linda J.
Review of General Psychology, Vol 10(3), Sep 2006, 229-240)
