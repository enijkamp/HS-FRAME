function [Gcos, Gsin, symbol] = GaborFilter2(scale, orient, shad)


% generate Gabor filter at fixed scale and orientation
% "Gcos" is the Gabor cos part
% "Gsin" is the Gabor sin part
% "symbol" is the bar for display
expand = 12; h = floor(scale * expand+.5);  
alpha = (pi * orient)/180;  % alpha is orientation
Gauss = zeros(h+h+1);  % Gaussian function
Gcos = zeros(h+h+1); Gsin = zeros(h+h+1); % Gabor cos and sine pair
symbol = zeros(h+h+1); 


for x0 = -h : h
    for y0 = -h : h
        if (x0^2+y0^2>h^2) 
            inCircle = 0;  % zeros for pixels outside the circle
        else
            inCircle = 1; 
        end
             
        x = (x0 * cos(alpha) + y0 * sin(alpha))/scale; 
        y = (y0 * cos(alpha) - x0 * sin(alpha))/scale;
                    
        g = exp(-(4*x^2+y^2)/100)/50/pi/scale^2; 
        Gauss(h+x0+1,h+y0+1) = g*inCircle; 
        Gcos(h+x0+1, h+y0+1) = g*cos(x)*inCircle;
        Gsin(h+x0+1, h+y0+1) = g*sin(x)*inCircle;
        
        %% make adaptive width of the bars 
        %symbol(h+x0+1, h+y0+1) = (abs(x)<3.4)*inCircle*shad; % make a bar
        %symbol(h+x0+1, h+y0+1) = (abs(x*2)<3.4)*inCircle*shad; % make a bar
        
        %% use uniform width of the bars
        symbol(h+x0+1, h+y0+1) = (abs(x*scale*1.7)<3.4)*inCircle*shad; % make a bar  % 1.4
        
    end
end

s = sum(Gauss(:)); sc = sum(Gcos(:)); r = sc/s; 
Gcos = Gcos - Gauss*r;   % mean is 0 by substracting DC component
Scos = sqrt(sum(sum(Gcos.^2))); Ssin = sqrt(sum(sum(Gsin.^2)));
Gcos = Gcos/Scos; Gsin = Gsin/Ssin; % l_2 norm is 1. 
%G = Gcos + sqrt(-1)*Gsin;  % sine and cosine pair







