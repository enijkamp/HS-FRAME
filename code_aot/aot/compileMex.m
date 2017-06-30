function [] = compileMex()

old_dir = cd('aot');

mex CsharedSketch.c;  % learning from aligned images
mex mexc_Sigmoid.cpp; % sigmoid transformation
mex Ccopy.c; % copy around detect location
mex ClocalNormalize.c
mex mexc_ComputeMAX1.cpp
mex mexc_ComputeSUM2.cpp
mex mexc_ComputeMAX2.cpp
mex mexc_ComputeSUM3.cpp
mex mexc_TemplateAffineTransform.cpp
mex mexc_CropInstance.cpp
mex mexc_FakeMAX2.cpp

cd(old_dir);

end

