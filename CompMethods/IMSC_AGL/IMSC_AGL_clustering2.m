function [Clu_result, indic] = IMSC_AGL_clustering2(X, truth, ind_folds, oracle, num_clusts)

if oracle
    numClust = length(unique(truth));
else
    numClust = num_clusts;
end

lambda1  = 0.1;
lambda2  = 1000;
lambda3  = 100;
numInst  = length(truth);
dim = numClust;

truthF = truth;  
numClust = length(unique(truthF));
num_view = length(X);
disp(num_view);

for iv = 1:num_view
    X1 = X{iv}';
    X1 = NormalizeFea(X1,1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];       
    Y{iv} = X1';            
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1;                                               
end
clear X X1 W1
X = Y;
clear Y      
for iv = 1:num_view
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 0;
    options.WeightMode = 'Cosine';
    Z1 = constructW(X{iv}',options);
    Z_ini{iv} = full(Z1);
    clear Z1;
end

for iv = 1:num_view
    invXX{iv} = inv(X{iv}'*X{iv}+2*eye(size(X{iv},2)));
end
F_ini = solveF(Z_ini,G,numClust);
U_ini = solveU(F_ini,numClust);

max_iter = 1;
miu = 0.01;
rho = 1.1;


[U] = IMSAGL(X,G,Z_ini,F_ini,invXX,numClust,lambda1,lambda2,lambda3,miu,rho,max_iter);
new_F = U;
norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
for i = 1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
new_F = new_F./norm_mat; 

indic    = kmeans(real(new_F),numClust);

result_CLU = ClusteringMeasure(truth, indic)*100;   

Clu_result.ACC = result_CLU(1);
Clu_result.NMI = result_CLU(2);
Clu_result.Purity = result_CLU(3);
% [Clu_result.ARi,~,~,~] = RandIndex(truth, indic);
[Clu_result.ARi] = 0;

end
