function [allS3SelectedRow, allS3SelectedCol, allS3SelectedOri]= rotateS3Template(S3SelectedRow,S3SelectedCol,S3SelectedOri,rotateShiftLimit,nOrient,numRotate,numPart)

allS3SelectedRow = zeros(numRotate,numPart,'single');
allS3SelectedCol = zeros(numRotate,numPart,'single');
allS3SelectedOri = zeros(numRotate,numPart,'single');


for (rot = -rotateShiftLimit : rotateShiftLimit)
    
        theta = -pi*rot/nOrient;
        sintheta = sin(theta);
        costheta = cos(theta);
            
        r = rot+rotateShiftLimit+1;  
        
        if (rot == 0)
            for (i=1: numPart)
                allS3SelectedRow(r,i) = S3SelectedRow(i);
                allS3SelectedCol(r,i) = S3SelectedCol(i);
                allS3SelectedOri(r,i) = S3SelectedOri(i);                      
            end
        else
            for (i=1: numPart)
                allS3SelectedRow(r,i) = S3SelectedRow(i)*costheta + S3SelectedCol(i)*sintheta;
                allS3SelectedCol(r,i) = -S3SelectedRow(i)*sintheta + S3SelectedCol(i)*costheta;          
                allS3SelectedOri(r,i) = S3SelectedOri(i) + rot;                       
            end       
        end       
end


     