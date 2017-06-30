Readme
======

Part I: Reference
-----------------

If you use the code please kindly cite the following paper:

[1] Jianwen Xie, Yifei Xu, Erik Nijkamp, Ying Nian Wu, Song-Chun Zhu. "Generative Hierarchical Learning of Sparse FRAME Models". CVPR 2017

The code is may be to reproduce the results shown in [1].


Part II: Dependencies
---------------------

The code requires the Matlab mex compiler to be configured for the compilation of C.

The package includes the external libraries liblinear 2.11-1 (multi-core variant) and piotrs toolbox 3.50.


Part III: How to run the code?
------------------------------

The code is configured to re-run the experiments published in [1]:

./config.m                 -  Contains parameters for 'AOT', see [1].
./main_1_codebook.m        -  Runs EM-algorithm to learn codebook. 
./main_2_features.m        -  Runs feature extraction step. 
./main_3_classication.m    -  Runs SVM classification step. 
./dataset                  -  Contains 'AnimalFace' dataset.
./aot/                     -  Contains code for AoT model.
./output/                  -  Contains the results.

Perform the following steps to run the experiments:

(1) config.m

- Set the configuration 'para.name' to 'cluster_5_toy', 'cluster_5', or 'cluster_11', e.g. para.name='cluster_5'.

- Set the number of cpu-cores 'para.noWorkers' to the number of CPU-threads, e.g. 'para.noWorkers = 8'.

(2) main_1_codebook.m 

- Run this file.

(3) main_2_features.m

- Run this file.

(4) main_3_classification.m

- Run this file.

