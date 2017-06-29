function logZRatio = computeLogZRatio(prevSamples,currSamples,filters,gradientF,step_width,n_iter)
% here we assume: 1. boundaries of gradientF is already set to 0;
% 2. prevSamples and currSamples are of the same sizes
%% copy items to GPU
gFilters = cell(size(filters));
for iFilter = 1:numel(filters)
    gFilters{iFilter}=gpuArray(filters{iFilter});
end
gI0 = gpuArray(prevSamples);
gI1 = gpuArray(currSamples);
[sx sy]=size(gradientF{1});
M = size(prevSamples,1)/sx;
N = size(prevSamples,2)/sy;
gDLambda= cell(size(gradientF));
for iFilter = 1:numel(gradientF)
    gDLambda{iFilter}=step_width*gpuArray(gradientF{iFilter});
end
%% compute the z ratio according to eqn 33.

logL0 = computeLogL(gI0,gFilters,sx,sy,M,N,gDLambda);
logL1 = computeLogL(gI1,gFilters,sx,sy,M,N,gDLambda);
logZRatio = mean([logL0,logL1]);
disp( [ 'LogZRatio: ' num2str(logZRatio,'%e') ]);
end

% L = q_1(I)/q_0(I)
function logLMat=computeLogL(I,gFilters,sx,sy,M,N,gDLambda)
    LMat =gpuArray.zeros(M,N);
   for iFilter = 1:numel(gFilters)
       rMap = abs(filter2(gFilters{iFilter},I));
       for m = 1:M
           for n = 1:N
               delta = rMap((m-1)*sx+1:m*sx,(n-1)*sy+1:n*sy).*gDLambda{iFilter};
               LMat(m,n)=LMat(m,n) + sum(delta(:));
           end
       end
   end
%   disp(LMat)
   LMat = reshape(gather(LMat),M*N,1);
   max_L = max(LMat);
   logLMat = max_L + log( sum( exp(LMat-max_L) ) )-log(M*N);

%   disp(logLMat)
end
