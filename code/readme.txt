Readme
======

Part I: Reference
-----------------

If you use the code please kindly cite the following paper:

[1] Jianwen Xie, Yifei Xu, Erik Nijkamp, Ying Nian Wu, Song-Chun Zhu. "Generative Hierarchical Learning of Sparse FRAME Models". CVPR 2017

The code is may be to reproduce the results shown in [1].


Part II: Dependencies
---------------------

The code requires the following Matlab toolboxes to be installed:

- Parallel Computing Toolbox
- Neural Network Toolbox

In addition, the Matlab mex compiler must be configured for the compilation of C.

The package includes the external libraries liblinear 2.11-1 (multi-core variant) and piotrs toolbox 3.50.


Part III: How to run the code?
------------------------------

The code is configured to re-run the experiments published in [1]:

./config.m                 -  Contains parameters for 'ours w/o parts' (baseline) and 'ours' (full), see [1].
./main_1_clustering.m      -  Runs clustering step. 
./main_2_classification.m  -  Runs classification step. 
./data/                    -  Contains 'AnimalFace' dataset.
./hs-frame/                -  Contains Code to train HS-FRAME model.
./output/                  -  Contains the results.

Perform the following steps to run the experiments:

(1) config.m

- Set the configuration 'para.name' to 'cluster_5_toy', 'cluster_5_full', 'cluster_5_base' or 'cluster_11_full', e.g. para.name='cluster_5_full'.

- Set the number of cpu-cores 'para.noWorkers' to the number of CPU-threads, e.g. 'para.noWorkers = 8'.

- ToDo

(2) main_1_clustering.m 

- ToDo

- Note, running the clustering step for a single category of the 'AnimalFace' dataset (i.e. config 'cluster_5_toy') will take about 72 CPU-hours (Intel Xeon E5-2676 v3). That is, the clustering task requires about 9 hours on a CPU with 8 threads / cores.


(3) main_2_classification.m

- ToDo

- Note, the seed for the pseudo-random number sequences used in the multi-core variant of liblinear is not fixed. Therefore, multiple executions of the code will result in results with slight variations.

