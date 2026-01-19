function [Clu_result, clust_lab]= OPIMC_clustering2(X, truth, ind_folds, bs, oracle, num_clusts)

if oracle
    numClust = length(unique(truth));
else
    numClust = num_clusts;
end

num_view = length(X);
numInst  = length(truth);

if size(X{1},2)~=numInst
    for iv = 1:num_view
        X{iv} = X{iv}';
    end
end
for iv = 1:length(X)
    X1 = X{iv};
    X1 = NormalizeFea(X1,0);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = 0 ;
    Y{iv} = X1; 
    W{iv} = diag(ind_folds(:,iv));                       
end

label = truth;
ind = ind_folds;
index = randperm(length(label));
for i = 1:num_view
%     X{i} = X{i}';           
    X{i} = X{i}(:,index);   
    W{i} = ind(index,i);    
end 
label = label(index);
block_size = bs;
option.label = label;
option.k       = numClust;
option.maxiter = 30;
option.tol     = 1e-6;
option.num_cluster = numClust;
option.pass = 1;
option.loss = 0;
option.alpha = 10;
[Clu_result, clust_lab] = OPIMC(X, W, option, block_size);
end
