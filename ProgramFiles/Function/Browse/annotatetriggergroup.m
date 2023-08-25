function annotatetriggergroup(msdata,scoredata,current_displaydataind,currentisomerind,...
    scoreintdata,theofragsto,currentCombiES,currentDecoyES)
[~,ind1] = sort(scoreintdata.colnames);
[~,ind2] = sort(scoreintdata.scoreoptions.analyzefragmode);
[~,ind3] = sort(ind1);
ind = ind2(ind3);
ms2tol = scoreintdata.scoreoptions.ms2tol(ind);
ms2tolunit = scoreintdata.scoreoptions.ms2tolunit(ind);
currentind = zeros(size(current_displaydataind));
numannofig = length(currentind);
for ii = 1:length(current_displaydataind)
    currentind(ii) = current_displaydataind{ii}(currentisomerind);
end
current_scoredata= scoredata(currentind(currentind > 0));
h = figure('Name',['CombiES = ',num2str(round(currentCombiES,4)),...
    '  CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
    '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
    '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein],...
    'Units','Pixels','NumberTitle','off',...
    'Position',[50 50 430*size(current_displaydataind,2) 720],'visible','off');
set(h,'unit','normalized','Tag','MainFigure'); % set figure as non-resizable and unit as normalized
% set uicontrol default value
set(h,'defaultuicontrolunits','normalized'); % set uicontrol default value unit as normalized
set(h,'defaultuicontrolfontsize',8);  % set uicontrol default value fontsize as 11
set(h,'defaultuicontrolfontname','Arial'); % set uicontrol default value font name as arial
set(h,'defaultuicontrolhorizontal','left'); % set uicontrol default value character horizontal as left

tabgp = uitabgroup(h,'OuterPosition',[0,0,1,1]);
tab1 = uitab(tabgp,'Title','Spectrum');
tab2 =  uitab(tabgp,'Title',['MS',char(178),' Detail']);
tab3 = uitab(tabgp,'Title',['MS',char(178),' Detail']);

figposfactor = 0.98/numannofig;
sgp = current_scoredata(1).SGP;
spectraplot = cell(1,numannofig);
errorplot = cell(1,numannofig);
textsummary = cell(1,numannofig);
legendplot = cell(1,numannofig);
foundpeaktable_handle = cell(1,numannofig);

spectra = msdata.spectra;
scannum = msdata.scannum;
charge = msdata.charge;

for ii = 1:numannofig
    if currentind(ii) == 0
        set(thisfigsto{ii},'visible','off');
    else
        
        
    end
    spectraplot{ii} = axes('parent',tab1,'Units','normalized','OuterPosition',...
        [figposfactor*(ii-1),0,figposfactor*.85,.45],...
        'TickDir','out','XLim',[0,2000]);
    errorplot{ii} = axes('parent',tab1,'Units','normalized','OuterPosition',...
        [figposfactor*(ii-1),.45,figposfactor*.85,.1],...
        'XLim',[0,2000]);
    textsummary{ii} = uitable('parent',tab1,'Units','normalized','OuterPosition',...
        [figposfactor*(ii-1)+.6/numannofig,.6,figposfactor*.4,.4],...
        'RowName',[]);
    foundpeaktable_handle{ii} = uitable('parent',tab2,'Units','normalized',...
        'OuterPosition',[figposfactor*(ii-1),.01,figposfactor*.98,.48]);
    thisscannum = current_scoredata(ii).Scan;
    tempspec = spectra{scannum == thisscannum};
    thischarge = charge(scannum == thisscannum);
    fragmethod = upper(current_scoredata(ii).Fragmode);
    denoisingoptions.fragmode = scoreintdata.scoreoptions.analyzefragmode;
    denoisingoptions.cutoffmed = scoreintdata.scoreoptions.cutoffmed;
    denoisingoptions.fracmax = scoreintdata.scoreoptions.fracmax;
    thisspec = spectradenoising(tempspec,fragmethod,thischarge,...
        glypepMW(sgp),ms2tol(ii),ms2tolunit{ii},denoisingoptions);
    spectrum = thisspec;
    theofrag = theofragsto{ii}{currentisomerind};
    cisoptions.maxlag = scoreintdata.scoreoptions.maxlag;
    cisoptions.selectpeak = scoreintdata.scoreoptions.selectpeak;
    matchresult = calcithscore(thisspec,[theofrag.mz],thischarge,...
        ms2tol(ii),ms2tolunit{ii},2,cisoptions);
    peakMatchIndex_merge=~cellfun(@isempty,matchresult.peakmatchindex);
    [isoIndex,~]=findIsotopePeaks(spectrum,ms2tol(ii),ms2tolunit{ii});
    for jj = 1:length(isoIndex)
        if any(peakMatchIndex_merge(isoIndex.mono(jj)))
            peakMatchIndex_merge(isoIndex(jj).iso,1)=1;
            peakMatchIndex_merge(isoIndex(jj).iso,2)=1;
            peakMatchIndex_merge(isoIndex(jj).iso,3)=1;
        else
            peakMatchIndex_merge(isoIndex(jj).iso,1)=0;
            peakMatchIndex_merge(isoIndex(jj).iso,2)=0;
            peakMatchIndex_merge(isoIndex(jj).iso,3)=0;
        end
    end
    width1=0.5;
    bar(spectraplot{ii},spectrum(sum(peakMatchIndex_merge,2) == 0,1),...
        spectrum(sum(peakMatchIndex_merge,2) == 0,2),...
        width1,'k','EdgeColor','k','FaceColor','k');
    % b/B-ion: blue
    bar(spectraplot{ii},spectrum(logical(peakMatchIndex_merge(:,3)),1),...
        spectrum(logical(peakMatchIndex_merge(:,3)),2),...
        width1,'b','EdgeColor','b','FaceColor','b');
    % y/Y-ion: red
    bar(spectraplot{ii},spectrum(logical(peakMatchIndex_merge(:,2)),1),...
        spectrum(logical(peakMatchIndex_merge(:,2)),2),...
        width1,'r','EdgeColor','r','FaceColor','r');
    % all others: green
    bar(spectraplot{ii},spectrum(logical(peakMatchIndex_merge(:,1)),1),...
        spectrum(logical(peakMatchIndex_merge(:,1)),2),...
        width1,'FaceColor',[30,106,21]/128,'EdgeColor',[30,106,21]/128);
    legendplot{ii} = axes('parent',tab1,'Units','normalized','Position',...
        [figposfactor*(ii-1)+.85/numannofig,0,.15*figposfactor,.35]);
    hold(legendplot{ii},'on');
    axis(legendplot{ii},'equal')
    set(legendplot{ii},'XLim',[0,5]);
    axis(legendplot{ii},'off')
    set(get(spectraplot{ii},'XLabel'),'String','mz');
    set(spectraplot{ii},'FontSize',12,'fontweight','bold');
    set(get(spectraplot{ii},'YLabel'),'String','Ion Current');
    ylim =get(spectraplot{ii},'ylim');
    set(spectraplot{ii},'ylim',[ylim(1) ylim(2)*1.16]);
    set(spectraplot{ii},'ClippingStyle','rectangle')
    foundPeaks = easymatch(spectrum,sgp,matchresult,theofrag);
end


end



























