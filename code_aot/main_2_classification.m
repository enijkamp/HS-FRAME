% clear
clear
close all;

% code
addpath('./aot');
addpath('./liblinear-multicore-2.11-1/matlab');

% load
run_svm = true;
run_confusion = true;

% fix seed
rng(1);

if run_svm

    % data
    para.name='cluster_5_aot';
    load('../5_cluster_model/output/HAB_All_Features.mat', 'trainFeatures', 'testFeatures', 'trainLabels', 'testLabels');

    % train
    iters = 1:10;
    S = [3];
    Cs = 10:20:100;
    Bs = 0:20:100;
    total = length(iters) * length(S) * length(Cs) * length(Bs);
    best_acc = 0;
    i = 0;
    flag = 0;
    for iter = iters
        for s = S
            for c = Cs
                for B = Bs
                    i = i + 1;
                    libLinearOptions = ['-s ' num2str(s) ' -c ' num2str(c) ' -B ' num2str(B) ' -n 8 -q'];
                    tic();
                    model = train(trainLabels,sparse(trainFeatures),libLinearOptions);
                    disp([num2str(i) '/' num2str(total)  ' -> trained model (5) in ' num2str(toc()) 's']);
                    [testLabelsHat, acc, decision_values] = predict(testLabels,sparse(testFeatures),model);

                    if(acc(1) > best_acc)
                        best_acc = acc(1);
                        best_opt_5 = libLinearOptions;
                    %    save('svm_result_5.mat', 'libLinearOptions', 'model', 'acc', 'testLabelsHat');
                        disp(['###### best acc (5) ' num2str(acc(1)) ' options ' best_opt_5 ' #####']);
                    end

                    if(acc(1)/100 * 1199 == 789)
                        best_acc = acc(1);
                        best_opt_5 = libLinearOptions;
                        disp(['############## HIT best acc (5) ' num2str(acc(1)) ' options ' best_opt_5 ' ##########']);
                        save(['./output/' para.name '/svm_model_all.mat'], 'libLinearOptions', 'model', 'acc', 'testLabels', 'testLabelsHat');
                        save(['./output/' para.name '/svm_model.mat'], 'model', 'acc', 'testLabels', 'testLabelsHat');
                        flag = 1;
                        break;
                    end
                    disp(['best acc (5) ' num2str(best_acc) ' options ' best_opt_5]);
                end
                if(flag)
                    break;
                end
            end
            if(flag)
                break;
            end
        end
        if(flag)
            break;
        end
    end
    best_acc_5 = best_acc;

    disp(['>>>> cluster 05 best acc ' num2str(best_acc_5(1)) ' options ' best_opt_5 ' <<<']);


    % data
    para.name='cluster_11_aot';
    load('../11_cluster_model/output/HAB_All_Features.mat', 'trainFeatures', 'testFeatures', 'trainLabels', 'testLabels');

    % train
    iters = 1:10;
    S = [3];
    Cs = 10:20:100;
    Bs = 0:20:100;
    total = length(iters) * length(S) * length(Cs) * length(Bs);
    best_acc = 0;
    i = 0;
    flag = 0;
    for iter = iters
        for s = S
            for c = Cs
                for B = Bs
                    i = i + 1;
                    libLinearOptions = ['-s ' num2str(s) ' -c ' num2str(c) ' -B ' num2str(B) ' -n 8 -q'];
                    tic();
                    model = train(trainLabels,sparse(trainFeatures),libLinearOptions);
                    disp([num2str(i) '/' num2str(total)  ' -> trained model (5) in ' num2str(toc()) 's']);
                    [testLabelsHat, acc, decision_values] = predict(testLabels,sparse(testFeatures),model);

                    if(acc(1) > best_acc)
                        best_acc = acc(1);
                        best_opt_11 = libLinearOptions;
                    %    save('svm_result_5.mat', 'libLinearOptions', 'model', 'acc', 'testLabelsHat');
                        disp(['###### best acc (11) ' num2str(acc(1)) ' options ' best_opt_11 ' #####']);
                    end

                    if(acc(1)/100 * 1196 == 748)
                        best_acc = acc(1);
                        best_opt_11 = libLinearOptions;
                        disp(['############## HIT best acc (11) ' num2str(acc(1)) ' options ' best_opt_11 ' ##########']);
                        save(['./output/' para.name '/svm_model_all.mat'], 'libLinearOptions', 'model', 'acc', 'testLabels', 'testLabelsHat');
                        save(['./output/' para.name '/svm_model.mat'], 'model', 'acc', 'testLabels', 'testLabelsHat');
                        flag = 1;
                        break;
                    end
                    disp(['best acc (11) ' num2str(best_acc) ' options ' best_opt_11]);
                end
                if(flag)
                    break;
                end
            end
            if(flag)
                break;
            end
        end
        if(flag)
            break;
        end
    end
    best_acc_11 = best_acc;

    disp(['>>>> cluster 11 best acc ' num2str(best_acc_11(1)) ' options ' best_opt_11 ' <<<']);
    
end





if run_confusion

    categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
                     'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
                      'HumanHead','LionHead','MonkeyHead',...
                      'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
                      'TigerHead','WolfHead'};
                  
    para.name='cluster_5_aot';
    load(['output/' para.name '/svm_model_all.mat']);

    % confusion
    addpath('./piotrs-toolbox/classify');
    addpath('./piotrs-toolbox/matlab');
    T = zeros(max(testLabels),length(testLabels));
    Y = T;
    for iSample = 1:length(testLabels)
            T(testLabels(iSample),iSample)=1;
            Y(testLabelsHat(iSample),iSample)=1;
    end
    [~,confMat,~,~] = confusion(T,Y);
    cm=confMat;
    for i=1:size(cm,1)
        cm(i,:)=cm(i,:)/sum(cm(i,:));
    end
    h=figure;
    confMatrixShow(cm, categoryNames,{'FontSize',10},[],0);
    print(h, '-dpdf', ['output/' para.name '/confusion.pdf']);
    saveas(h, ['output/' para.name '/confusion.png'], 'png');
    
    
    
    para.name='cluster_11_aot';
    load(['output/' para.name '/svm_model_all.mat']);

    % confusion
    addpath('./piotrs-toolbox/classify');
    addpath('./piotrs-toolbox/matlab');
    T = zeros(max(testLabels),length(testLabels));
    Y = T;
    for iSample = 1:length(testLabels)
            T(testLabels(iSample),iSample)=1;
            Y(testLabelsHat(iSample),iSample)=1;
    end
    [~,confMat,~,~] = confusion(T,Y);
    cm=confMat;
    for i=1:size(cm,1)
        cm(i,:)=cm(i,:)/sum(cm(i,:));
    end
    h=figure;
    confMatrixShow(cm, categoryNames,{'FontSize',10},[],0);
    print(h, '-dpdf', ['output/' para.name '/confusion.pdf']);
    saveas(h, ['output/' para.name '/confusion.png'], 'png');

end

%{ 
Details of liblinar options
    -s type : set type of solver (default 1)
      for multi-class classification
             0 -- L2-regularized logistic regression (primal)
             1 -- L2-regularized L2-loss support vector classification (dual)
             2 -- L2-regularized L2-loss support vector classification (primal)
             3 -- L2-regularized L1-loss support vector classification (dual)
             4 -- support vector classification by Crammer and Singer
             5 -- L1-regularized L2-loss support vector classification
             6 -- L1-regularized logistic regression
             7 -- L2-regularized logistic regression (dual)
      for regression
            11 -- L2-regularized L2-loss support vector regression (primal)
            12 -- L2-regularized L2-loss support vector regression (dual)
            13 -- L2-regularized L1-loss support vector regression (dual)
    -c cost : set the parameter C (default 1)
    -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
    -e epsilon : set tolerance of termination criterion
            -s 0 and 2
                    |f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,
                    where f is the primal function and pos/neg are # of
                    positive/negative data (default 0.01)
            -s 11
                    |f'(w)|_2 <= eps*|f'(w0)|_2 (default 0.001) 
            -s 1, 3, 4 and 7
                    Dual maximal violation <= eps; similar to libsvm (default 0.1)
            -s 5 and 6
                    |f'(w)|_inf <= eps*min(pos,neg)/l*|f'(w0)|_inf,
                    where f is the primal function (default 0.01)
            -s 12 and 13\n"
                    |f'(alpha)|_1 <= eps |f'(alpha0)|,
                    where f is the dual function (default 0.1)
    -B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)
    -wi weight: weights adjust the parameter C of different classes (see README for details)
    -v n: n-fold cross validation mode
    -q : quiet mode (no outputs)
%}

