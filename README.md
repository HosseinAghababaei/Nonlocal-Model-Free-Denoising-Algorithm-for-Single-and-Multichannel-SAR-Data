In this repository you can find the matlab code associated to the IEEE TGRS paper 
[Nonlocal Model-Free Denoising Algorithm for Single- and Multichannel SAR Data](https://ieeexplore.ieee.org/document/9611274)).

This is meant as a unified framework for SAR denosining: single channel, multi-channel and interferometric phase
The code has meant for reaerch purpose
If you use it, please cite as the following:

H. Aghababaei, G. Ferraioli, S. Vitale, R. Zamani, G. Schirinzi and V. Pascazio, "*Nonlocal Model-Free Denoising Algorithm for Single- and Multichannel SAR Data,*" 
in IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-15, 2022, Art no. 5217315, doi: 10.1109/TGRS.2021.3127109.

# Usage 
You can find:
1. the matlab code of the method
```
Stc_McSAR.m
```
2. A free available ESAR testing image
```
ESAR.mat
```
3. Auxiliar codes
```
mask_window.m
weight_fcn.m
Pauli_C.m
```
4. A demo for running the code on the available data:
```
demo.m
```
In order to filter your specified the code you should provide as input:
- **Cin:** the covariance matrix of the image
- **Patch_size:** the local window size (typical values are 3x3)
- **win_size:** the search window size (typical values are 21x21 or 23x23)
- **Look:** is the number of look used for estimating the covariance matrix

for more details, please refer to the paper implementation details or contact us.

# Team members
Hossein Aghababaei (contact person, h.aghababaei@utwente.nl)
Giampaolo Ferraioli (giampaolo.ferraioli@uniparthenope.it);
Sergio Vitale    (sergio.vitale@uniparthenope.it);
Roghayeh Zamani
Gilda Schirinzi (gilda.schirinzi@uniparthenope.it)
Vito Pascazio (vito.pascazio@uniparthenope.it)
 
# License
Copyright (c) 2022 Department of Earth Observation Science, Faculty of Geo-Information Science and Earth Observation, University of Twente
aand Dipartimento di Ingegneria and Dipartimento di Scienze e Tecnologie of Università degli Studi di Napoli "Parthenope".

All rights reserved. This work should only be used for nonprofit purposes.

By downloading and/or using any of these files, you implicitly agree to all the
terms of the license, as specified in the document LICENSE.txt
(included in this directory)

