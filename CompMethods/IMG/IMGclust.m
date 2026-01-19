function [Ux Uy P2 P1 P3 F P R nmi avgent AR C] = IMGclust(X2,Y2,X1,Y3,numClust,truth,option,oracle,num_clusts)

if (min(truth)==0)
    truth = truth + 1;
end

if oracle
    numClust = length(unique(truth));
else
    numClust = num_clusts;
end

option.option = numClust;
option.truth = truth;
[Ux,Uy,P2,P1,P3,W] = IMG(X2,Y2,X1,Y3,option);

% fprintf('running spectral clustering...\n');
kmeans_avg_iter = 10;
for i=1: kmeans_avg_iter

    C = clu_ncut(W,numClust);
    C = C';
    
    %%
    [A nmii(i) avgenti(i)] = compute_nmi(truth,C);
    [Fi(i),Pi(i),Ri(i)] = compute_f(truth,C);
    [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
end
F = mean(Fi);
P = mean(Pi);
R = mean(Ri);
nmi = mean(nmii);
avgent = mean(avgenti);
AR = mean(ARi);

fprintf('nmi: %f(%f)\n', nmi, std(nmii));

