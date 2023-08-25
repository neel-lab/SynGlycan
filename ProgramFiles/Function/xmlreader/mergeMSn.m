function [msdata] = mergeMSn(msdata,simscoretol) % merge all MSn spectra that in each cluster that are similar
j = waitbar(0,'','Name','Performing Merging...');
allscannum = msdata.scannum;
sim_score = msdata.simScore;
consensus = msdata.spectra;
cluster_idx = cell(length(msdata.spectra),1);
for ii = 1 : length(sim_score)
    scan_idx = find(allscannum == sim_score{ii,1});
    sim_score1 = sim_score{ii,2};
    if ~isempty(sim_score1)
        for jj = 1:length(sim_score1(:,1))
            if sim_score1(jj,2) >= simscoretol
                sim_scan_idx = allscannum == sim_score1(jj,1);
                combine_spec = [consensus{scan_idx};msdata.spectra{sim_scan_idx}]; % add spectra in a cluster
                [~,idx] = sort(combine_spec(:,1)); % sort combined spectrum based on m/z
                consensus{scan_idx} = combine_spec(idx,:);
                consensus{sim_score1(jj,1)} = []; % removing spectra for sim_score1(jj,1)
                cluster_idx{scan_idx} = [cluster_idx{scan_idx},sim_score1(jj,1)];
            end
        end
    end
    waitbar(ii/length(sim_score),j,sprintf('%d/%d',ii,length(sim_score)));
end
msdata.spectra = consensus;
msdata.clusterID = cluster_idx;
close(j)
end