function logZRatio = computeLogZRatio_sparse(prevSamples,currSamples,template, filters,gradientF_selected,step_width,nRow, nCol, sx, sy)
% here we assume: 1. boundaries of gradientF is already set to 0;
% 2. prevSamples and currSamples are of the same sizes
%% copy items to GPU
gI0 = prevSamples;
gI1 = currSamples;

% gDLambda= cell(size(gradientF));
% for iFilter = 1:numel(gradientF)
%     gDLambda{iFilter}=step_width*gpuArray(gradientF{iFilter});
% end

DLambda = step_width*gradientF_selected;

%% compute the z ratio according to eqn 33.

logL0 = computeLogL(gI0,template,filters,nRow, nCol,DLambda, sx, sy);
logL1 = computeLogL(gI1,template,filters,nRow, nCol,DLambda, sx, sy);
logZRatio = mean([logL0,logL1]);
disp( [ 'LogZRatio: ' num2str(logZRatio,'%e') ]);
end

% L = q_1(I)/q_0(I)
function logLMat=computeLogL(I,template, filters,nRow, nCol, DLambda, sx, sy)
     
   
   
%    rSample = sparseFeatureExtraction(single(I),filters, pTemplate,true);
%    
%    LMat = zeros(nRow, nCol);
%    for i = 1:nRow*nCol
%        delta = rSample(1+(i-1)*numSelected : numSelected +(i-1)*numSelected).*DLambda;
%        LMat(i) = sum(delta(:));
%    end
   

   LMat = zeros(nRow, nCol);
   for i = 1:nRow
      for j=1:nCol
           singleChain_I = I(  1+(i-1)*sx : sx+(i-1)*sx, 1+(j-1)*sx : sy+(j-1)*sy); 
           rSample = sparseFeatureExtraction( single(singleChain_I ),filters, template, true);  
           delta = rSample.*DLambda;
           LMat(i) = sum(delta(:));
      end
   end


   LMat = reshape(LMat,nRow*nCol,1);
   max_L = max(LMat);
   logLMat = max_L + log( sum( exp(LMat-max_L) ) )-log(nRow*nCol);   % log ( (exp(x_1)+...+exp(x_M))/M ),  where LMat=x_1, ..., x_M. 

%   disp(logLMat)
end
