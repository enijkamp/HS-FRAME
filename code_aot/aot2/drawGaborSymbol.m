function im = drawGaborSymbol(im, allsymbol, row, col, orientationIndex, nGaborOri, scaleIndex, intensity)
% drawGaborSymbol - Display one Gabor symbol on image canvas.

h = floor( (size(allsymbol{(scaleIndex-1)*nGaborOri + orientationIndex}, 1)-1)/2 );  % half size of Gabor
for r = row-h:row+h
    if r < 1 || r > size(im,1)
        continue;
    end
    for c = col-h:col+h
        if c < 1 || c > size(im,2)
            continue;
        end
        % disp(sprintf('ori=%d,xx=%d,yy=%d',orientationIndex,r-row+h+1,c-col+h+1));
        val = intensity * allsymbol{(scaleIndex-1)*nGaborOri + orientationIndex}(r-row+h+1,c-col+h+1);
        if val > im(r,c)
            im(r,c) = val;
        end
    end
end
