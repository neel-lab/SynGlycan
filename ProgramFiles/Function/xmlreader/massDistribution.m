function theomassdist = massDistribution(scannumind,mz,chg,allretime,allmslevel,allspectra,defaultmonoisodist,chnos)
dist_range = defaultmonoisodist{2,1}(1,1)-defaultmonoisodist{1,1}(1,1);
mass = mz * chg - 1.0078246*chg;
tempretime = allretime(scannumind);
temprange = (abs(allretime-tempretime) <= .25) & allmslevel == 1;
spectrain=allspectra(temprange);
spectrasize = cellfun(@size,spectrain,num2cell(ones(size(spectrain)))); % calculate size of each of these spectra
spectra = zeros(sum(spectrasize),2); % create empty element with size equal to the sum of the elements in all spectra
ind = 1;
for i = 1:length(spectrasize)        % merge all input spectra in a single variable.
    spectra(ind:ind + spectrasize(i) - 1,:) = double(spectrain{i});
    ind = ind + spectrasize(i);
end
spectrumout = spectra(spectra(:,2) > 0,:);
spectrumsectionind = ((spectrumout(:,1) - mz) <= 5/chg) & ((spectrumout(:,1) - mz) >= -3/chg);
spectrumsection = spectrumout(spectrumsectionind,:);  % consider chg >= 2, 3 Da off monoiso max
idx = floor(mass);
if idx == 0 % if the precursor m/z is 0
    error('No precursor m/z and/or charge value present in msdata.')
elseif idx <= defaultmonoisodist{end,1}(1,1) % if the precursor m/z is within the theoretical distribution
    idx_alt = round(idx/dist_range); % find index of mass based on distribution range
    if idx_alt == 0 % if rounding results in 0
        idx_alt = 1; % set value to 1
    elseif idx_alt > length(defaultmonoisodist) % if rounding results in value greater than length of theoretical distribution
        idx_alt = length(defaultmonoisodist); % set to max
    end
    theoisodist = defaultmonoisodist{idx_alt,1}; % select theoretical distribution
else
    % [0.024321959339170,0.476559048284606,0.001461698046926,0.012941719677186,0] GLYCAN CHNOS VALUES
    % [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292] GLYCOPEPTIDE CHNOS VALUES
    if ~isempty(chnos)
        avggp = chnos;
    else
        avggp = [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292]; % [C H N O S]
    end
    [theoisodist,~,~] = isotopicdist(mass*avggp,'ShowPlot',0);
end
mzfix = (mass - floor(mass)) / chg;
theoisodist(:,1) = theoisodist(:,1)/chg + 1.0078246 + mzfix;
theoisostep = mean(diff(theoisodist(:,1)));
expandtheoiso = flip(theoisodist(1,1)-theoisostep:-theoisostep:min(spectrumsection(:,1)));
expandtheoiso = [expandtheoiso(:),zeros(length(expandtheoiso),1)];
expandtheoiso = [expandtheoiso;theoisodist];
etisotail = theoisodist(end,1)+theoisostep:theoisostep:max(spectrumsection(:,1))+1;
expandtheoiso = [expandtheoiso;[etisotail(:),zeros(length(etisotail),1)]];
expandtheoiso = unique(expandtheoiso,'rows');
% Collecting all peaks in expandtheoiso about 10%
totIC = sum(expandtheoiso(:,2));
theomassdist = [expandtheoiso(:,1),(expandtheoiso(:,2)/totIC)];
intidx = any(theomassdist(:,2) < 0.1,2);
theomassdist(intidx,:) = [];
end