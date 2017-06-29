function [purity, entropy] = evaluateClustering(inPath, nameGT, probs)

% purity = \sum_c p(c) \max_{l\in Gt} p(l|c)
purity = 0;

% entropy =  \sum_c p(c) \sum_{l\in Gt} p(l|c) log {1\over p(l|c)}
entropy = 0;

try
    fid = fopen([inPath '/' nameGT], 'r');
catch
   warning(' can not find groundtruth file\n'); 
   return;
end

tmp = textscan(fid, '%s %d'); % imgfilename gtLabel

gt.imgFileNames = tmp{1};
gt.label        = tmp{2};

% imgFiles used to do EM
imgFiles = dir([inPath '/*.jpg']);
numImg   = length(imgFiles);
label    = zeros(numImg, 1);
for i = 1 : numImg
    idx = strcmp(imgFiles(i).name, gt.imgFileNames);
    if sum(idx(:))~=1
        error(' can not find groundtruth for %s', imgFiles(i).name);
    end        
    label(i) = gt.label(idx);
end

k = size(probs, 2);
clusterImgIdx = cell(k,1);
for i = 1:size(probs,1)
    [~, ind] = max(probs(i,:));    
    clusterImgIdx{ind}=[clusterImgIdx{ind}, i];
end

% compute
numImg = length(label);
labelVal = unique(label);
numLVal = length(labelVal);

for c = 1 : k
    curClusterNum = length(clusterImgIdx{c});
    curClusterGt  = label(clusterImgIdx{c});
    
    pc = curClusterNum / numImg;
    
    plc = zeros(1, numLVal);
    entropylc = zeros(1, numLVal);
    for i = 1 : numLVal
        plc(i) = sum(curClusterGt==labelVal(i)) / curClusterNum;
        if plc(i)>0
            entropylc(i) = plc(i) * log(1/(plc(i)));
        end
    end    
    
    purity = purity + pc * max(plc(:));
    
    entropy = entropy + pc * sum(entropylc(:));
end
