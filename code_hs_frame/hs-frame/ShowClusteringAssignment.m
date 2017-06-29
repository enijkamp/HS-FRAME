    %%%%% generating html
    fid = fopen(fullfile(resultPath, ['result' num2str(it) '.html']), 'wt');
    fprintf(fid, 'Original images<br><br>\n');
    for c=1:numCluster
       for i=1:size(clusters(c).imageIndex,2)
           id=clusters(c).imageIndex(i);
           fprintf(fid, '%s\n', ['<IMG SRC="' fullfile('./img', files(id).name) '" height=' num2str(sx)  ' width=' num2str(sy) '>']); 
       end
       fprintf(fid, '%s\n', ['<br>']);
       
       fprintf(fid, '%s\n', ['<hr>']);
    end

    fprintf(fid, '%s\n', ['<br>']);
    fprintf(fid, '%s\n', ['<hr>']);
    
    fprintf(fid, 'Morphed images<br><br>\n');
    for c=1:numCluster
       for i=1:size(clusters(c).imageIndex,2)
           id=clusters(c).imageIndex(i);
           fprintf(fid, '%s\n', ['<IMG SRC="' fullfile( ['./iteration' num2str(it)], 'morphedCropped', ['morphed-cluster-' num2str(c) '-img-' num2str(id,'%04d') '.png']) '" height=' num2str(sx)  ' width=' num2str(sy) '>']); 
       end
       fprintf(fid, '%s\n', ['<br>']);
       
       fprintf(fid, '%s\n', ['<hr>']);
    end
    
    if it~=0
    fprintf(fid, '%s\n', ['<br>']);
    fprintf(fid, '%s\n', ['<hr>']);
    
    fprintf(fid, 'images with bounding box<br><br>\n');
    for c=1:numCluster
       for i=1:size(clusters(c).imageIndex,2)
           id=clusters(c).imageIndex(i);
           fprintf(fid, '%s\n', ['<IMG SRC="' fullfile( ['./iteration' num2str(it)], 'boundingBox', ['boundingBox-cluster-' num2str(c) '-img-' num2str(id,'%04d') '.png']) '" height=' num2str(sx)  ' width=' num2str(sy) '>']); 
       end
       fprintf(fid, '%s\n', ['<br>']);
       
       fprintf(fid, '%s\n', ['<hr>']);
    end
    end
    fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
