%
% Learning mixture of And-Or Template Model from Non-aligned Images
%

clear

if (exist('output','dir'))
    delete('output/*.*'); 
else
    mkdir('output');
end

if ~exist('working','dir')
    mkdir('working');
end
for c = 1:5
    mkdir(['output/' num2str(c)]);
end
Compile;

categoryss = cell(12,1);
categoryss{1} = '2_bull_cow';
categoryss{2} = '2_cup_teapot';
categoryss{3} = '2_plane_helicopter';
categoryss{4} = '3_camel_elephant_deer';
categoryss{5} = '3_clock';
categoryss{6} = '3_swan_eagle_seagull';
categoryss{7} = '4_eye_mouth_ear_nose';
categoryss{8} = '4_flowers';
categoryss{9} = '4_pc';
categoryss{10} = '5_deer_cat_wolf_tiger_lion';
categoryss{11} = '5_musical_intrument';
categoryss{12} = '5_tiger_leopard_zebra_dog_cow';
NumClusters = [2 2 2 3 3 3 4 4 4 5 5 5];
total_purity = zeros(12, 5);
total_entropy = zeros(12,5);


for cate = 11  % specify the the task ID
    for seed = 1:4
        clearvars -except categoryss NumClusters seed cate total_purity total_entropy;
        delete('output/*.*'); 
        loadInitialization = true;
        NumCluster = NumClusters(cate);
        categorys = categoryss{cate};
        SetParameter;
        numIteration = 15;
        SUM1MAX1;
        EM_Clustering; 
        total_purity(cate,seed) = purity(numIteration);
        total_entropy(cate, seed) = entropy(numIteration);
    end
end

total_purity
total_entropy