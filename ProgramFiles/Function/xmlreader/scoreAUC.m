function sim_score = scoreAUC(scan_clusters,allspectra,MS2tol,MS2tolUnit) % scoring scan clusters
g = waitbar(0,'','Name','Performing Similarity Scoring...');
if (strcmpi(MS2tolUnit,'Da')) % for Da
    MS2tolC = MS2tol;  % mass tolerance in Da
elseif ~(strcmpi(MS2tolUnit,{'ppm','Da'}))
    error('MATLAB:GLYCOPAT:ERRORUNIT','INCORRECT UNIT FOR MS2 MASS TOL');
end
[~,c] = size(scan_clusters);
if c > 4
    aa = unique(scan_clusters(:,5));
    if aa(1) == 0
        aa = aa(2:end);
    end
    sim_score = cell(length(aa),2);
    for mm = 1:length(aa)
        scan_idx = (scan_clusters(:,5) == aa(mm));
        ri = scan_clusters(scan_idx,:);
        lvl = 1;
        [rows,~] = size(ri);
        for rr = 1:rows
            row = nonzeros(ri(rr,:));
            if length(row) > 4
                if strcmpi(MS2tolUnit,'ppm')  % for ppm
                    MS2tolC = MS2tol/1e6*row(1);  % mass tolerance in ppm
                end 
                scans = row(4:end);
                lscn = length(scans);
                scannum = row(5);
                sim_score{mm,1} = scannum;
                raw_spec1 = ms2centroid(allspectra{scannum}); % centroid spectrum
                if isempty(raw_spec1)
                    continue
                end
                [~,idx1] = sort(raw_spec1(:,2),'descend');
                spec1 = raw_spec1(idx1,:);
                m1 = padarray(spec1(:,1),30,'post');
                i1 = padarray(spec1(:,2),30,'post');
                for nn = 1:lscn
                    if nn == 2
                        continue
                    end
                    raw_spec2 = ms2centroid(allspectra{scans(nn)});
                    if isempty(raw_spec2)
                        continue
                    end
                    % Similarity Score Calculation
                    [~,idx2] = sort(raw_spec2(:,2),'descend');
                    spec2 = raw_spec2(idx2,:);
                    m2 = padarray(spec2(:,1),30,'post');
                    i2 = padarray(spec2(:,2),30,'post');
                    for k = 1:30
                        hp(k) = i1(k); % part numerator
                        hp_m = m1(k);
                        mat = abs(m2-hp_m) <= MS2tolC;
                        if any(mat)
                            loc = find(mat);
                            if length(loc) > 1
                                [~,loc] = min(abs(m2-hp_m));
                            end
                            hpi(k) = i2(loc); % part numerator
                        else
                            hpi(k) = 0; % part numerator
                        end
                        numer(k) = hp(k)*hpi(k); % total numerator
                        hp_tot(k) = (hp(k))^2; % part denominator
                        hq(k) = (i2(k))^2; % part denominator
                    end
                    sim_score{mm,2}(lvl,1) = scans(nn);
                    sim_score{mm,2}(lvl,2) = sum(numer)/sqrt(sum(hp_tot)*sum(hq)); % similarity score calculation
                    lvl = lvl+1;
                end
            end      
        end
        waitbar(mm/length(aa),g,sprintf('%d/%d',mm,length(aa)));
    end
else
    sim_score = [];
end
close(g)
end