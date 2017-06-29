function [rModel_selected, pSample]=multiChainHMC_sparse(template, filters, pCurrSample,stepsize,L, nRuns, nRow, nCol)
% function to sample from multiple parallel chains, with gaussian prior
% nRow and nCol defines the size of the chain grid.
% This is a CPU version of sparse multiple chain HMC sampling 
% we assume the pCurrSample is already a large canvas containing a grid of nRow by nCol sampled images

%% 0. prepare the enlarged template
pTemplate = template.enlargedTemplate;   % p for paralle
pTemplate.selectedLambdas = repmat(template.selectedLambdas,1,nRow*nCol);
numSelected = template.numSelected;


%% 1. correct the bounary condition, and run HMC
pSample = pCurrSample;
for iRun= 1:nRuns
    pSample = multiHMC(pTemplate,filters,pSample,stepsize,L,nRow,nCol, numSelected);
end
%% 2 compute the rModel_selected

[pRSample] = sparseFeatureExtraction(single(pSample),filters, pTemplate,true);
rModel_selected = mean(reshape(pRSample,numSelected,[]),2);
rModel_selected = rModel_selected';

end % end of main function


function [q]=multiHMC(pTemplate, filters, pCurrSample, stepsize, L, nRow, nCol, numSelected)
% The main routine for sampling using HMC.
%
%% The traditional leapFrog steps

q=pCurrSample;
sigma2 = .1; % sigma^2

p = randn(size(q),'single');
p = p*sqrt(sigma2);

current_p=p;

p=p- (stepsize/2) * grad_U_sparse(q, filters, pTemplate);


for i=1:L
  q=q+ stepsize * p/sigma2;
  if(i~=L) 
      p=p-stepsize * grad_U_sparse(q, filters, pTemplate); 
  end
end

p=p-stepsize/2* grad_U_sparse(q, filters, pTemplate);


%% compute the accept probability for each chain separately 
[currentUMap_I, currentUMap_main]= multiUMap_sparse(pCurrSample, filters, pTemplate, numSelected, nRow, nCol);

[proposedUMap_I, proposedUMap_main]= multiUMap_sparse(q, filters, pTemplate, numSelected, nRow, nCol);


current_p = current_p.^2; % reuse the variable, to save GPU memory
p = p.^2;
[sx sy]=size(pCurrSample);
sx = sx/nRow; sy = sy/nCol;  % now sx and sy is the size of each small image
keys = rand(nRow,nCol);
for iRow = 1:nRow
    for iCol = 1:nCol
        sx0 = (iRow-1)*sx+1;
        sx1 = iRow*sx;
        sy0 = (iCol-1)*sy+1;
        sy1 = iCol*sy;
        
        currentU = sum(sum(currentUMap_I(sx0:sx1,sy0:sy1)))+currentUMap_main(iRow,iCol);
        currentK = sum(sum(current_p(sx0:sx1,sy0:sy1)))/2/sigma2;
        
        proposedU = sum(sum(proposedUMap_I(sx0:sx1,sy0:sy1)))+proposedUMap_main(iRow,iCol);
        proposedK = sum(sum(p(sx0:sx1,sy0:sy1)))/2/sigma2;
        
        accProb = exp(currentU-proposedU+currentK-proposedK);
        %fprintf('accept ratio for chain (%d,%d) is %f, ',iRow,iCol,accProb);
        if keys(iRow,iCol)<accProb
            %fprintf('accepted!\n');
        else
            %fprintf('rejected!\n');
            q(sx0:sx1,sy0:sy1)=pCurrSample(sx0:sx1,sy0:sy1);
        end
    end %iRow
end % iCol

end % gpuMultiHMC
        

      
function [pEnergyUMap_I, pEnergyMap_main] = multiUMap_sparse(pCurrImage, filters, pTemplate, numSelected, M,N)
   
 
   [rSample] = sparseFeatureExtraction(single(pCurrImage),filters, pTemplate,true);
   
   pEnergyMap_main=reshape(-pTemplate.selectedLambdas(1:numSelected)*reshape(rSample,numSelected,[]),M,N);
   
   pEnergyUMap_I = pCurrImage.^2/2;
   
      
end


