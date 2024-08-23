# A deep learning-based crystal plasticity finite element model

This is the code for paper "A deep learning-based crystal plasticity finite element model"

## Installation Requirements
The basic requirement for using the files is a Python 3.6.3 environment with PyTorch 2.3.0
## Source Files

- MLtrain.py is the code for model training.
- Cp-ML-PowerLaw-func.py is the main code file for the whole CP-AI framework.
- Others are related functions of the CP-AI framework.

## Running the code
If you want to predict a stress-strain curve using this framework, please change the parameters in input file py.15, py_init.15, texture.15 and texture_init.15 (two init files include parameters for initial CP steps, the other two files include parameter for ML models). Then, run `CP-ML-PowerLaw-func.py`. The output is stress-strain data saved in file SS.txt.

## Developer Team
The code was developed by Yuwei Mao from the [CUCIS](http://cucis.ece.northwestern.edu/index.html) group at the Electrical and Computer Engineering Department at Northwestern University.

## Disclaimer
The research code shared in this repository is shared without any support or guarantee on its quality. However, please do raise an issue if you find anything wrong and I will try my best to address it.

email: yuweimao2019@u.northwestern.edu

Copyright (C) 2023, Northwestern University.

See COPYRIGHT notice in top-level directory.

## Funding Support
This work is supported in part by the following grants: National Institute of Standards and Technology (NIST) award 70NANB19H005; Department of Energy (DOE) award DE-SC0021399; Na- tional Science Foundation (NSF) awards CMMI-2053929, OAC-2331329; and Northwestern Center for Nanocombinatorics.
