function [Crr, Cri, Cir, Cii] = CorrFilter(Fr, Fi)

% C = corr(F)
%---------------------------------------------------------
% Compute correlation or inhibition coefficients so that 
% after the synthesis coefficient of a base is updated, 
% the analysis coefficients of all the overlapping bases 
% can be updated without computing the remainder image. 
%---------------------------------------------------------

% F: the bank of filters
% C: the correlation matrix

N = size(Fr, 2);    % N: the number of filters

for (i = 1:N) 
   hi = floor((size(Fr{i}, 1)-1)/2+.5);   % half-size of filter i
   for (j = 1:N)
      hj = floor((size(Fr{j}, 1)-1)/2.+.5);    % half-size of filter j
      Ir = zeros(2*(hi+hj)+1);   
      Ir((hj+1):(hj+2*hi+1), (hj+1):(hj+2*hi+1)) = Fr{i}; 
      Ii = zeros(2*(hi+hj)+1);   
      Ii((hj+1):(hj+2*hi+1), (hj+1):(hj+2*hi+1)) = Fi{i}; 
        
      
      Crr{i, j} = filter2(Fr{j}, Ir, 'same');
      Cri{i, j} = filter2(Fi{j}, Ir, 'same');
      Cir{i, j} = filter2(Fr{j}, Ii, 'same');
      Cii{i, j} = filter2(Fi{j}, Ii, 'same');
      
   end
end



