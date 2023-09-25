# CalibrateWALOP

## Introduction
This is a software suite to create a model of the WALOP instrument,
as well as a map of measured standard deivation of the inversion of the model
for different polarization states.

## Prerequisites
- Python 3.9
- astropy >= 5.1
- matplotlib >= 3.6.0
- numpy >= 1.23.3
- pandas >= 1.5.0
- pip >= 22.2.2
- pipenv >= 2022.9.8
- plotly >= 5.11.0
- scipy >= 1.9.1
- sympy >= 1.11.1

## Usage
### Preparation of files
In this directory place the following tree:
```
{This folder}
│
└───Config 1
│   │   transmission_field1_pol1.txt
│   │   transmission_field1_pol2.txt
│   │   transmission_field1_pol3.txt
│   │   ...
│   │   transmission_field1_pol49.txt
│   │   transmission_field2_pol1.txt
│   │   transmission_field2_pol2.txt
│   │   transmission_field2_pol3.txt
│   │   ...
│   │   transmission_field2_pol49.txt
│   │   transmission_field3_pol1.txt
│   │   transmission_field3_pol2.txt
│   │   transmission_field3_pol3.txt
│   │   ...
│   │   ...
│   │   transmission_field12_pol1.txt
│   
└───Config 2
│   │   transmission_field1_pol1.txt
│   │   ...
│   │   ...
│   │   transmission_field12_pol1.txt
│   
└───Config 2
│   │   transmission_field1_pol1.txt
│   │   ...
│   │   ...
│   │   transmission_field12_pol1.txt
│   
└───Config 2
│   │   transmission_field1_pol1.txt
│   │   ...
│   │   ...
│   │   transmission_field12_pol1.txt
│   
│   ...
```
Each 'Config' folder must contain exactly 
12(field ensembles)x49(polarization states)x12(fields/ensemble) 
accoring to the methodology in Maharana et. al. 2022, Kypriotakis et. al.
2023 (in prep).

A sample such file is shown below:
```agsl
Polarization Transmission Data

File : C:\Users\johny\OFFLINE\EngDrwV3.5.1.zmx
Title: WALOP-North
Date : 20/4/2022
Configuration 1 of 4

X-Field     :  1.0000
Y-Field     :  0.0000
X-Phase     :    0.00
Y-Phase     :    0.00

Grid Size : 32 x 32

Aperture, Fresnel, coating, vignetting, and internal transmittance effects are considered.

Field Pos : -0.2500, -0.2500 (deg)
  Transmission at    0.5400:     0.096654699
  Transmission at    0.5550:     0.094061271
  Transmission at    0.5800:     0.090250045
  Transmission at    0.6200:     0.092560855
  Transmission at    0.6600:     0.099341866
  Transmission at    0.6850:     0.097649242
  Transmission at    0.7000:     0.092031494
  Total Transmission       :     0.094770944

Field Pos : -0.2500, -0.2000 (deg)
  Transmission at    0.5400:     0.101592624
  Transmission at    0.5550:     0.097104901
  Transmission at    0.5800:     0.091935702
  Transmission at    0.6200:     0.095772629
  Transmission at    0.6600:     0.100970087
  Transmission at    0.6850:     0.099786485
  Transmission at    0.7000:     0.094441957
  Total Transmission       :     0.097117560

Field Pos : -0.2500, -0.1600 (deg)
  Transmission at    0.5400:     0.104581171
  Transmission at    0.5550:     0.101600013
  Transmission at    0.5800:     0.094805969
  Transmission at    0.6200:     0.097488384
  Transmission at    0.6600:     0.103244632
  Transmission at    0.6850:     0.101665980
  Transmission at    0.7000:     0.096469159
  Total Transmission       :     0.099764040

Field Pos : -0.2500, -0.1100 (deg)
  Transmission at    0.5400:     0.105452061
  Transmission at    0.5550:     0.102403457
  Transmission at    0.5800:     0.097600535
  Transmission at    0.6200:     0.097830453
  Transmission at    0.6600:     0.104360103
  Transmission at    0.6850:     0.102362041
  Transmission at    0.7000:     0.097367697
  Total Transmission       :     0.100913304

Field Pos : -0.2500, -0.0700 (deg)
  Transmission at    0.5400:     0.105874696
  Transmission at    0.5550:     0.102795707
  Transmission at    0.5800:     0.097899723
  Transmission at    0.6200:     0.100561598
  Transmission at    0.6600:     0.104882214
  Transmission at    0.6850:     0.102673521
  Transmission at    0.7000:     0.097792605
  Total Transmission       :     0.101762836

Field Pos : -0.2500, -0.0200 (deg)
  Transmission at    0.5400:     0.105745364
  Transmission at    0.5550:     0.102662173
  Transmission at    0.5800:     0.097731331
  Transmission at    0.6200:     0.100919840
  Transmission at    0.6600:     0.102830811
  Transmission at    0.6850:     0.102442743
  Transmission at    0.7000:     0.097659660
  Total Transmission       :     0.101318914

Field Pos : -0.2500, 0.0200 (deg)
  Transmission at    0.5400:     0.105649515
  Transmission at    0.5550:     0.102564166
  Transmission at    0.5800:     0.097635367
  Transmission at    0.6200:     0.100804411
  Transmission at    0.6600:     0.102710752
  Transmission at    0.6850:     0.102330487
  Transmission at    0.7000:     0.097565696
  Total Transmission       :     0.101210624

Field Pos : -0.2500, 0.0700 (deg)
  Transmission at    0.5400:     0.105555602
  Transmission at    0.5550:     0.102465557
  Transmission at    0.5800:     0.097571766
  Transmission at    0.6200:     0.100174242
  Transmission at    0.6600:     0.104470598
  Transmission at    0.6850:     0.102279237
  Transmission at    0.7000:     0.097465100
  Total Transmission       :     0.101392750

Field Pos : -0.2500, 0.1100 (deg)
  Transmission at    0.5400:     0.104975567
  Transmission at    0.5550:     0.101902980
  Transmission at    0.5800:     0.097093453
  Transmission at    0.6200:     0.097260731
  Transmission at    0.6600:     0.103706585
  Transmission at    0.6850:     0.101759046
  Transmission at    0.7000:     0.096874984
  Total Transmission       :     0.100346872

Field Pos : -0.2500, 0.1600 (deg)
  Transmission at    0.5400:     0.103940244
  Transmission at    0.5550:     0.100907732
  Transmission at    0.5800:     0.094127930
  Transmission at    0.6200:     0.096640176
  Transmission at    0.6600:     0.102301242
  Transmission at    0.6850:     0.100832060
  Transmission at    0.7000:     0.095808899
  Total Transmission       :     0.098965464

Field Pos : -0.2500, 0.2000 (deg)
  Transmission at    0.5400:     0.100874813
  Transmission at    0.5550:     0.096342258
  Transmission at    0.5800:     0.091116435
  Transmission at    0.6200:     0.094711592
  Transmission at    0.6600:     0.099822544
  Transmission at    0.6850:     0.098811506
  Transmission at    0.7000:     0.093696853
  Total Transmission       :     0.096165349

Field Pos : -0.2500, 0.2500 (deg)
  Transmission at    0.5400:     0.095910820
  Transmission at    0.5550:     0.093178326
  Transmission at    0.5800:     0.089229919
  Transmission at    0.6200:     0.091230066
  Transmission at    0.6600:     0.097963925
  Transmission at    0.6850:     0.096541571
  Transmission at    0.7000:     0.091232644
  Total Transmission       :     0.093628534


Chief Ray Transmission Surface By Surface:

Field Pos : -0.2500, -0.2500 (deg)
Wavelength 1: 0.540 �m

Surf    Tot. Tran    Rel. Tran
 Chief ray does not trace or is vignetted.

Field Pos : -0.2500, -0.2500 (deg)
Wavelength 2: 0.555 �m

Surf    Tot. Tran    Rel. Tran
 Chief ray does not trace or is vignetted.

...
...
...
```

### Run the code
If you have completed the preparation, run the files run_inv_1089_nonoise.py 
and run_inv_1089_noise.py.

If all went well you must have the files inv_1089_nonoise.csv and 
inv_1089_noise.csv generated.

Then run the file plotter_part.py, changing line 12:
```python
with open('inv_1089_nonoise.csv') as csv_file:
```

according to your needs.