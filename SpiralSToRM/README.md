# SpiralSToRM-Iterative/SpiralSToRM-Navigator

main file: SpiralStorm_iterative/SpiralStorm_navigator

Input Kpace Data: Number of Frequency points x number of channels x number of spirals 

Spiral SToRM Iterative/Spiral SToRM Navigator 
Reference paper:

A.H. Ahmed, Y. Mohsin, R. Zhao, Y. Yang, M. Salerno, P. Nagpal, M. Jacob, 
Free-breathing and ungated cardiac cine using navigator-less spiral SToRM
(submitted to MRM)

Link:https://arxiv.org/pdf/1901.05542.pdf

In the above paper, we have developed an algorithm for a free-breathing and ungated cardiac MRI scheme using an iterative kernel low-
rank algorithm and self-gated spiral sequence.

Main benefits:

The  iterative  SToRM  algorithm  facilitates  the  extension  of  manifold  regularization  to
navigator-less spiral acquisitions, thus improving sampling efficiency. The proposed scheme eliminates
the need for breath-holding and ECG gating, while facilitating the imaging of the cardiac function in
different respiratory phases.

Dependencies
Using gpuNUFFT provided by:
gpuNUFFT - GPU Regridding of Radial 3-D MRI data
 - Andreas Schwarzl - andreas.schwarzl@student.tugraz.at
 - Florian Knoll - florian.knoll@nyumc.org

Dataset

We have released two dataset(approx~only 5 sec of data) used in this paper. You can download the full dataset from the below link:

Download Links : https://drive.google.com/open?id=1oscN89LaaHcDlT5O23y3GyMmOlQd-Fqu
                https://drive.google.com/open?id=1Z8o7p-OFDFatt3DVEW6wQs0YA2RolW99
This dataset consist of Cardiac cine MRI: k space data: Frequency encoding pts x channels x no.of spirals

The code is provided to support reproducible research. 
If the code is giving any error or some files are missing then you may open an issue 
or directly email me at abdul-ahmed@uiowa.edu or hasib_bhati@yahoo.com