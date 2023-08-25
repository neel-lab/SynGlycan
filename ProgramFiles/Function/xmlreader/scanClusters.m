function scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz)
h = waitbar(0,'','Name','Generating Spectra Clusters...');
for ww = 1:length(allscannum)
    if allmslvl(ww)>1 % for all MSn scans
        retime = allretime(ww);
        precursorMz = allprecursormz(ww);
        origparse = scan_clusters(:,1);
        difference = abs(origparse - precursorMz) <= MS1tolC; % look if precursor m/z is already in matrix within mass tolerance
        if any(difference) % if precursor m/z is already in matrix
            kk = find(difference); % row of precursor m/z if already exists
            for jj = 1:length(kk) % for each row that precursor m/z is in
                if allscannum(ww) == scan_clusters(kk(jj),4) % if scan number is already present, skip
                    continue
                elseif (retime >= scan_clusters(kk(jj),2)) && (retime <= scan_clusters(kk(jj),3)) % if precursor m/z is within retention range
                    [r,c] = size(scan_clusters);
                    massLocation = zeros(r,1);
                    scan_clusters = [scan_clusters, massLocation]; % add column of zeros to matrix
                    d = c + 1;
                    for ll = 1:d
                        if scan_clusters(kk(jj),ll) == 0
                            scan_clusters(kk(jj),ll) = allscannum(ww); % replace first zero value in row with scan number
                        break
                        end
                    end
                    idx = (scan_clusters(:,4) == allscannum(ww));
                    scan_clusters(idx,:) = zeros;
                end
            end          
        end
    end
    waitbar(ww/length(allscannum),h,sprintf('%d/%d',ww,length(allscannum)));
end
scan_clusters(~any(scan_clusters,2),:) = []; % remove all rows with only zero values
scan_clusters(:,~any(scan_clusters,1)) = []; % remove all columns with only zero values
close(h)
end