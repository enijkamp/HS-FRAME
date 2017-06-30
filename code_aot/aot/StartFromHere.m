%
% LEARNING HIERARCHICAL ACTIVE BASIS MODEL FROM NONALIGNED IMAGES
%


clear;

if (exist('output','dir'))
    delete('output/*.*'); 
else
    mkdir('output');
end

if ~exist('working','dir')
    mkdir('working');
end

categorys = '5_tiger_leopard_zebra_dog_cow';
loadInitialization = true;
NumCluster = 5;
for c = 1:NumCluster
    mkdir(['output/' num2str(c)]);
end

SetParameter;
%Compile;

SUM1MAX1;

EM_Clustering; 