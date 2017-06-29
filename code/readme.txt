Readme
======

Part I: Reference
-----------------

If you use the code please kindly cite the following paper:

[1] Jianwen Xie, Yifei Xu, Erik Nijkamp, Ying Nian Wu, Song-Chun Zhu. "SGenerative Hierarchical Learning of Sparse FRAME Models". CVPR 2017

The code is may be to reproduce the results shown in [1].


Part II:  How to use the code?
------------------------------



requires:
- parallel computing toolbox
- neural network toolbox
- mex
- seed liblinear-multicore

config
- set number of cpu cores
- note running a single category of the 'animalface' dataset will require about X cpu-hours, that is with e.g. 8 cores it takes ... cluster_5_toy config

(1)
select model parameters in config.m

(2)
main_1

(3)
main_2








********************************
Part II:  How to use the code?
********************************

Please use the paremeters described in the paper of [1]

(0) Setup and compile: run setup.m 

(1) Experiment 1, 2, and 3.  (synthesis experiments)

Start from main.m.   The file for parameters setting is config_STGConvNet.m

%%%   Model type:
%%%   ST_3:         3 convolutional layer (layer by layer) (Exp 1 in the paper)
%%%   FC_S_3_large: 1 convolutional layer + 2 spatial fully connected layer (end to end)  (Exp 2 in the paper 224x224)
%%%   FC_S_3:       1 convolutional layer + 2 spatial fully connected layer (end to end)  (Exp 2 in the paper)
%%%   FC_ST_2:      1 convolutional layer + 1 spatial-temporal fully connected layer (end to end)  (Exp 3 in the paper)
%%%
%%%   e.g.,  
%%%   (1) Exp 1: category='sea';   type='ST_3'; 
%%%   (2) Exp 2: category='water_fountain';   type='FC_S_3_large'; 
%%%   (3) Exp 3: category='fire_pot';   type='FC_S_3'; 
%%%   (4) Exp 4: category='tiger';   type='FC_ST_2'; 

(2) large size of training data.  (synthesis experiments with mini-batch version)

start from main_largeScale.m.    The file for parameters setting is config_STGConvNet_largeScale.m

(3) Experiment 4 (recovery experiment)

start from main_recovery.m.    The file for parameters setting is config_STGConvNet_recovery.m

(4) Experiment 5 (Background inpainting)

start from main_background_inpainting.m.    The file for parameters setting is config_STGConvNet_background_inpainting.m

