function logZRatio = computeLogZRatio_sparse(prevSamples,currSamples,template, filters,gradientF_selected,step_width,nRow, nCol)
% here we assume: 1. boundaries of gradientF is already set to 0;
% 2. prevSamples and currSamples are of the same sizes
%% copy items to GPU
pTemplate = template.enlargedTemplate;
numSelected = template.numSelected;
gI0 = prevSamples;
gI1 = currSamples;

% gDLambda= cell(size(gradientF));
% for iFilter = 1:numel(gradientF)
%     gDLambda{iFilter}=step_width*gpuArray(gradientF{iFilter});
% end

DLambda = step_width*gradientF_selected;

%% compute the z ratio according to eqn 33.

logL0 = computeLogL(gI0,pTemplate,filters,nRow, nCol,DLambda,numSelected);
logL1 = computeLogL(gI1,pTemplate,filters,nRow, nCol,DLambda,numSelected);
logZRatio = mean([logL0,logL1]);
disp( [ 'LogZRatio: ' num2str(logZRatio,'%e') ]);
end

% L = q_1(I)/q_0(I)
function logLMat=computeLogL(I,pTemplate, filters,nRow, nCol, DLambda,numSelected)
     
   
   
   rSample = sparseFeatureExtraction(single(I),filters, pTemplate,true);
   
   LMat = zeros(nRow, nCol);
   for i = 1:nRow*nCol
       delta = rSample(1+(i-1)*numSelected : numSelected +(i-1)*numSelected).*DLambda;
       LMat(i) = sum(delta(:));
   end
   
%    for iFilter = 1:numel(filters)
%        rMap = abs(filter2(filters{iFilter},I));
%        for m = 1:M
%            for n = 1:N
%                delta = rMap((m-1)*sx+1:m*sx,(n-1)*sy+1:n*sy).*gDLambda{iFilter};
%                LMat(m,n)=LMat(m,n) + sum(delta(:));
%            end
%        end
%    end
%   disp(LMat)
   LMat = reshape(LMat,nRow*nCol,1);
   max_L = max(LMat);
   logLMat = max_L + log( sum( exp(LMat-max_L) ) )-log(nRow*nCol);

%   disp(logLMat)
end
