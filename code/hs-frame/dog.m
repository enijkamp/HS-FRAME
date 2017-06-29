function [B, B0, C] = dog(sigma, shad)

c = floor(sigma*2 + .5);        % the half size
siz = c+c-1; 
B = zeros(siz);        % B: the base matrix 
C = zeros(siz); 
B0 = zeros(siz); 
r = .9; s = 1;          % the sizes of center and sorround 

for i = 1:siz
  for j = 1:siz
    x = (i-c)/sigma;
    y = (j-c)/sigma;    % scale x and y
    B(i, j) = G(x, r)*G(y, r) - G(x, s)*G(y, s);
    ra = sqrt(x*x+y*y); 
    if ((ra<=.6)&&(ra>=.4)) 
         C(i, j) = shad;       % a circle 
    end
  end
end

B = B - mean(B(:)); 
tmp = B(:); 
B = B/sqrt(sum(tmp.*tmp));  % normalize to have norm 1



