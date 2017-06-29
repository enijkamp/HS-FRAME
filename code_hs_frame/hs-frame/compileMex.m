function [] = compileMex()

old_dir = cd('hs-frame');

mex Ccopy.c;
mex CgetMAX1.c;
mex CgetMAX1_DoG.c;
mex filterTiling.c;
mex sparseFeatureExtraction.c;
mex ComputeCorrBB.c;
mex I_plus_cB.c;
mex gibbs_v3.c;
mex sparseFiltering.c;
mex mexc_updateCorrelationTable.c;
mex drawSketch.c;
mex sparseFRAME_SUM2_part.c;
mex mexc_sparseFRAME_MAX2.cpp;
mex MatchingPursuitV3.cpp;
mex grad_U_sparse.c;
mex mexc_ComputeSUM3_logZ_partial.cpp;
mex mexc_TemplateAffineTransform.cpp;

cd(old_dir);

end

