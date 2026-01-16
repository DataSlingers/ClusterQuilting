
clear all;
clc

% Change folder name to where the data are stored.
% Assumes that the folder is contains subfolders of different sim settings,
% which contain further subfolders of individual simulation iterations. 
%% NB: If the comparison methods code is in a different directory, 
%% the addpath lines below need to be changed as well.
foldname = "./";
ggg = dir(foldname);

for ii = 3:length(ggg)
    rr = getfield(ggg(ii), 'name');
    rrr = strcat(foldname, rr);
    ww = dir(rrr);
    disp(rr);
    for jj = 3:length(ww)
        
        %Load data
        qq = getfield(ww(jj), 'name');
        if contains(qq, '.csv')
            continue
        end
        if startsWith(qq, "._");
            continue
        end
        disp(jj);
        sobs = readtable(strcat(rrr, "/", qq, "/sobs.csv"));
        stimes = readtable(strcat(rrr, "/", qq, "/stimes.csv"));
        testdat = table2array(readtable(strcat(rrr, "/", qq, "/masked_dat.csv")));
        testomega = table2array(readtable(strcat(rrr, "/", qq, "/omega.csv")));
        clustass = table2array(readtable(strcat(rrr, "/", qq, "/clustass.csv")));
        pom = table2array(readtable(strcat(rrr, "/", qq, "/panel_obs_mat.csv")));
        ptm = table2array(readtable(strcat(rrr, "/", qq, "/panel_times_mat.csv")));
        
        numClust = length(unique(clustass));
        [bs, ~] = size(sobs);
        [~, num_view] = size(pom);
        [nobs, ~] = size(testdat);
        
        
        for vv = 1:num_view
            ghj = find(ptm(:, vv) == 1);
            testdict{vv} = testdat(:, ghj);
            testdict_t{vv} = testdat(:, ghj)';
        end
        
        %%%%%%%%%%
        %%%%%%%%%%
        %IMSC-AGL
        disp("IMSC-AGL")
        addpath('./code_dir/CompMethods/IMSC_AGL');
        addpath('./code_dir/CompMethods/clustering_metrics');
        addpath('./code_dir/CompMethods/utils');

        
        [imscagl_res, imscagl_clust] = IMSC_AGL_clustering2(testdict_t, clustass, pom);
        writematrix(imscagl_clust,strcat(rrr, "/", qq, "/IMSCAGL_res.csv"));
        
        rmpath('./code_dir/CompMethods/IMSC_AGL');
        rmpath('./code_dir/CompMethods/clustering_metrics');
        rmpath('./code_dir/CompMethods/utils');

        
        %%%%%%%%%%
        %%%%%%%%%%
        %IMG
        disp("IMG")
        addpath('./code_dir/CompMethods/IMG');
        addpath('./code_dir/CompMethods/IMG/measure');
        addpath('./code_dir/CompMethods/IMG/misc');
        for v3 = 2:num_view
            paired = find(pom(:, 1) == 1 & pom(:, v3) == 1);
            spared = find(pom(:, 1) == 1 | pom(:, v3) == 1);
            if ~isempty(paired)
                singleInstView1 = find(pom(:, 1) == 1 & pom(:, v3) == 0);
                singleInstView2 = find(pom(:, 1) == 0 & pom(:, v3) == 1);
                xpaired=testdict{1}(paired,:);
                ypaired=testdict{v3}(paired,:);
                xsingle=testdict{1}(singleInstView1,:);
                ysingle=testdict{v3}(singleInstView2,:);
                option.latentdim=numClust;
                option.lamda=1e-2;
                option.beta=1;
                option.gamma=1e2;
                truthF = clustass(spared);
                try
                    [U1 U2 P2 P1 P3 F P R nmi avgent AR img_clust] = IMGclust(xpaired,ypaired,xsingle,ysingle,numClust,truthF,option);
                catch
                    img_clust = ones(nobs);
                end
                break
            end        
        end
        
        writematrix(img_clust,strcat(rrr, "/", qq, "/IMG_res.csv"));
        clear option
        rmpath('./code_dir/CompMethods/IMG');
        rmpath('./code_dir/CompMethods/IMG/measure');
        rmpath('./code_dir/CompMethods/IMG/misc');
        
        %%%%%%%%%%
        %%%%%%%%%%
        
        %DAIMC
        disp("DAIMC")
        addpath('./code_dir/CompMethods/DAIMC');
        addpath('./code_dir/CompMethods/clustering_metrics');
        addpath('./code_dir/CompMethods/utils');
        
        [daimc_res, daimc_clust] = DAIMC_clustering2(testdict, clustass, pom);
        writematrix(daimc_clust,strcat(rrr, "/", qq, "/DAIMC_res.csv"));
        
        rmpath('./code_dir/CompMethods/DAIMC');
        rmpath('./code_dir/CompMethods/clustering_metrics');
        rmpath('./code_dir/CompMethods/utils');
        
        %%%%%%%%%%
        %%%%%%%%%%
        %OPIMC
        disp("OPIMC")
        addpath('./code_dir/CompMethods/OPIMC');
        addpath('./code_dir/CompMethods/clustering_metrics');
        addpath('./code_dir/CompMethods/utils');
        
        [opimc_res, opimc_clust] = OPIMC_clustering2(testdict, clustass, pom, bs);
        writematrix(opimc_clust,strcat(rrr, "/", qq, "/OPIMC_res.csv"));
        
        rmpath('./code_dir/CompMethods/OPIMC');
        rmpath('./code_dir/CompMethods/clustering_metrics');
        rmpath('./code_dir/CompMethods/utils');
        
        %%%%%%%%%%
        %%%%%%%%%%
        %Nuclear Norm
        disp("NN")
        addpath('/Volumes/T7 Shield/CQ');

        [n2, p2] = size(testdat);
        try
            [Lt, Sigmat,outputt] = nucmin_dat(testomega, testdat(testomega), 1/max(eig(corr(testdat))),100, 0.01, n2, p2);
        catch
            Sigmat = zeros(n2, p2);
        end
        writematrix(Sigmat,strcat(rrr, "/", qq, "/nn_impute.csv"));

        rmpath('/Volumes/T7 Shield/CQ');
    end
end