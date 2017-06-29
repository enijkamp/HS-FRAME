function [Corr] = CorrFilterFrame(filters)

% C = corr(F)
%---------------------------------------------------------
% Compute correlation or inhibition coefficients so that 
% after the synthesis coefficient of a base is updated, 
% the analysis coefficients of all the overlapping bases 
% can be updated without computing the remainder image. 
% This is a version for using Frame version of filters bank
%---------------------------------------------------------

% F: the bank of filters
% C: the correlation matrix

N = size(filters, 2);    % N: the number of filters

for (i = 1:N) 
   hi = floor((size(filters{i}, 1)-1)/2+.5);   % half-size of filter i (filter i is considered as being filtered)
   for (j = 1:N)
      hj = floor((size(filters{j}, 1)-1)/2.+.5);    % half-size of filter j (filter j is considered as a filter)
      I = zeros(2*(hi+hj)+1);   
      I((hj+1):(hj+2*hi+1), (hj+1):(hj+2*hi+1)) = filters{i}; 
     
      Corr{i, j} = filter2(filters{j}, I, 'same');
       
   end
end