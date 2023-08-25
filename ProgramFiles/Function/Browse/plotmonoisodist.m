function plotmonoisodist(msdata,current_scoredata,MS1dist,MS1disttext,tmipwindow)
% FUNCTION: function discription
% 
% 
% Syntax:
% 
% 
% Input:
% 
% 
% Output:
% 
% 
% Note:
% 
% 
% Example:
% 
% 
% Children function: 
% 
% 
% See Also:
% 

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

defaultisodist = msdata.defaultisodist;
if isempty(MS1dist)
    MS1dist = axes;
end
cla(MS1dist);
hold(MS1dist,'on');
allretime = msdata.retime;
allspectra = msdata.spectra;
allmslvl = msdata.mslvl;
allprecmz = msdata.precursormz;
allscannum = msdata.scannum;
allprecscannum = cell2mat(msdata.precursorScanNum);

resultscannum = current_scoredata(1).Scan;
resultcharge = current_scoredata(1).Charge;
resultscanind = allscannum == resultscannum;
resultprecscanind = allscannum == (allprecscannum(allscannum == resultscannum));
parentretime = allretime(resultprecscanind);
precursorMz = allprecmz(resultscanind);
originalprecursorMz = msdata.allprecursormz{resultscanind};
candimz = glypepMW(current_scoredata(1).SGP)/resultcharge + 1.007825032;

retimeind = allretime >= (parentretime - tmipwindow) & allretime  <= (parentretime + tmipwindow);
spectraformonoiso = allspectra(retimeind & allmslvl == 1);
spectrasize = cellfun(@size,spectraformonoiso,...
    num2cell(ones(size(spectraformonoiso)))); % calculate size of each of these spectra
monospec = zeros(sum(spectrasize),2);
% create empty element with size equal to the sum of the elements in all spectra
ind = 1;
for ii = 1:length(spectrasize)        % merge all input spectra in a single variable.
    monospec(ind:ind + spectrasize(ii) - 1,:) = double(spectraformonoiso{ii});
    ind = ind + spectrasize(ii);
end
spectrumout = monospec(monospec(:,2) > 0,:);
spectrumout = spectrumout((spectrumout(:,1) <= precursorMz + 2.5) &...
    (spectrumout(:,1) >= precursorMz - 2.5),:);
spectrumout = sortrows(spectrumout,2);
spectrumout = flipud(spectrumout);
[~,ia,~] = unique(spectrumout(:,1));
spectrumout = spectrumout(ia,:);
bar(MS1dist,spectrumout(:,1),spectrumout(:,2),'EdgeColor','k','FaceColor','k','barwidth',1);
pause(0.01);
MS1dist.XLabel.String = 'Isotope distribution  -  m/z';
MS1dist.YLabel.String = 'Intensity';
tempmaxint = max(spectrumout(:,2))*1.05;
if isempty(tempmaxint)
    tempmaxint = 1;
end
precmass = (precursorMz-1.007825032)*resultcharge;
thisdefaultisodist = defaultisodist{floor(precmass)};
thisdefaultisodist(:,1) = thisdefaultisodist(:,1)/resultcharge + 1.007825032;
thisdefaultisodist(:,1) = thisdefaultisodist(:,1) + candimz - thisdefaultisodist(1,1);
candipeakmzs = candimz + [-0.0125,0.0125];
candipeakint = spectrumout(spectrumout(:,1) >= candipeakmzs(1) & spectrumout(:,1) <= candipeakmzs(2),2);
if isempty(candipeakint)
    for ii = 2:size(thisdefaultisodist,1)
        candipeakmzs = candimz + 1.007825032/resultcharge*(ii - 1) + [-0.0125,0.0125];
        candipeakint = spectrumout(spectrumout(:,1) >= candipeakmzs(1) & spectrumout(:,1) <= candipeakmzs(2),2);
        if ~isempty(candipeakint)
            candipeakint = max(candipeakint)/thisdefaultisodist(ii,2)*thisdefaultisodist(1,2);
            break
        end
    end
end
thisdefaultisodist(:,2) = thisdefaultisodist(:,2) * (max(candipeakint))/(thisdefaultisodist(1,2));
plot(MS1dist,thisdefaultisodist(:,1),thisdefaultisodist(:,2),'-b');
scatter(MS1dist,candimz,max(candipeakint),'xk','SizeData',50);

plot(MS1dist,[candimz,candimz],[0,tempmaxint],'--g');
text(MS1dist,candimz,tempmaxint*.9,'Candi.','color','g')
plot(MS1dist,[originalprecursorMz(1),originalprecursorMz(1)],[0,tempmaxint],'r');
text(MS1dist,originalprecursorMz(1),tempmaxint*.8,'Prec.(Inst.)','color','r')
plot(MS1dist,[precursorMz,precursorMz],[0,tempmaxint],'b');
text(MS1dist,precursorMz,tempmaxint*.7,'Prec.(Avg.)','color','b');
auxisomarker = [fliplr(precursorMz:-1.007825032/resultcharge:precursorMz - 2.5),...
    precursorMz + 1.007825032/resultcharge:1.007825032/resultcharge:precursorMz + 2.5];
scatter(MS1dist,auxisomarker,zeros(size(auxisomarker)),'*m','SizeData',25);
set(MS1dist,'XLim',ceil(candimz) + [-2 2]);
if ~isempty(MS1disttext)
    if tmipwindow == 0
        set(MS1disttext,'String','Parent MS1 Scan');
    else
        set(MS1disttext,'String',['+/- ',num2str(tmipwindow),' min(s)']);
    end
end
end