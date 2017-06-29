function enlargedTemplate = enlargingTemplate(nTileRow, nTileCol, sx, sy, template)

    enlargedTemplate.selectedRow=[];
    enlargedTemplate.selectedCol=[];
    enlargedTemplate.selectedFilter=[];
    
   
   for i=1: nTileCol    % do not change the i and j order, which will affect the function of multiChainHMC_sparse.m
       for j=1:nTileRow
           for k=1:template.numSelected
           
           newRow=template.selectedRow(k) + (j-1) * sx;
           newCol=template.selectedCol(k) + (i-1) * sy;
           newFilter=template.selectedFilter(k);
           
           enlargedTemplate.selectedRow=[enlargedTemplate.selectedRow newRow];
           enlargedTemplate.selectedCol=[enlargedTemplate.selectedCol newCol];
           enlargedTemplate.selectedFilter=[enlargedTemplate.selectedFilter newFilter];
          
           end
       end
   end