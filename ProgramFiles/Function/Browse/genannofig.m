function genannofig(scannum,spectrum,score,MS2tol,MS2tolUnit,msdata,matchresult,...
    theofrag,scoretablefilename)
% GENANNOFIG: generate the 3 tab annotation figure from matching result
%
% Syntax:
% genannofig(scannum,spectrum,score,MS2tol,MS2tolUnit,msdata,matchresult,...
%     theofrag,scoretablefilename)
%
% Input:
% scannum: scan number of the spectrum
% spectrum: the experimental spectrum
% score: score returned by MATCH1BY1
% MS2tol: MS2 matching tolerence value.
% MS2tolUnit: MS2 matching tolerence unit.
% msdata: all experimental data.
% matchresult: result returned by CALCITHSCORE, must contain field
%     "peakmatchindex".
% theofrag: theoretical fragments of candidate (glyco)peptide.
% scoretablefilename: (if desired) where to save the details of matched
%     peaks.
%
% Output:
% The 3 tab figure.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% EASYMATCH
%
% See Also:
% GETDGFRAG  DRAWGLYCAN

sgp = score.SGP;
foundPeaks = easymatch(spectrum,sgp,matchresult,theofrag);
if isempty(foundPeaks)
    warndlg('No peak matched.');
    return
else
    textcolor = 'kbrmgyc';
    spacing = .01 * (max(spectrum(:,1)) - min(spectrum(:,1)));
    figdisplay = 1;
    tmipwindow = 0.25; % 0.5 min time window for monoisotopic pks. disp.
    [p,g,m] = breakGlyPep(sgp);
    allretime = msdata.retime;
    alltotIonCurrent = msdata.totIonCurrent;
    allmslvl = msdata.mslvl;
    allscannum = msdata.scannum;
    allspectra = msdata.spectra;
    scannumind = allscannum == scannum;
    precursorMz = msdata.precursormz(scannumind);
    precursorscannum = msdata.precursorScanNum{scannumind};
    originalprecursorMz = msdata.allprecursormz{scannumind};
    ms2retime = allretime(scannumind);
    parentretime = allretime(allscannum == precursorscannum);
    ismslvl1 = allmslvl == 1;
    
    %% Version 1. Plot Precursor for the entire experiment
    % thisTICdata = zeros(sum(ismslvl1),2);
    % ms1sn = allscannum(ismslvl1);
    % for i = 1:length(ms1sn)
    %     tempspec = allspectra{allscannum == ms1sn(i)};
    %     tempretime = allretime(allscannum == ms1sn(i));
    %     precpkind = abs(tempspec(:,1) - precursorMz) / precursorMz * 1e6 <= 10;
    %     thisTICdata(i,:) = [tempretime,sum(tempspec(precpkind,2))];
    % end
    
    %% BELOW IS A VERSION THAT ONLY PLOTS 20 MIN PRECURSOR TIC, IS FASTER
    retime10min = abs(allretime - parentretime) <= 5;
    thisTICdata = zeros(sum(ismslvl1 & retime10min),2);
    ms1sn = allscannum(ismslvl1 & retime10min);
    selectedspectra = allspectra(ismslvl1 & retime10min);
    selectedretime = allretime(ismslvl1 & retime10min);
    for i = 1:length(ms1sn)
        tempspec = selectedspectra{i};
        tempretime = selectedretime(i);
        precpkind = abs(tempspec(:,1) - precursorMz) / precursorMz * 1e6 <= 10;
        thisTICdata(i,:) = [tempretime,sum(tempspec(precpkind,2))];
    end
    
    plottpnum = 20;
    figdispname  = ['Glycoproteomics Analysis of ','Single Spectrum (scan number= ',...
        strtrim(int2str(scannum)),')'];
    %%set up figure properties
    h = figure('Name',figdispname,'Units','Pixels','NumberTitle','off',...
        'Position',[50 50 1280 720],'visible','off');
    
    set(h,'unit','normalized','Tag','MainFigure'); % set figure as non-resizable and unit as normalized
    % set uicontrol default value
    set(h,'defaultuicontrolunits','normalized'); % set uicontrol default value unit as normalized
    set(h,'defaultuicontrolfontsize',8);  % set uicontrol default value fontsize as 11
    set(h,'defaultuicontrolfontname','Arial'); % set uicontrol default value font name as arial
    set(h,'defaultuicontrolhorizontal','left'); % set uicontrol default value character horizontal as left
    
    tabgp = uitabgroup(h,'OuterPosition',[0,0,1,1]);
    tab1 = uitab(tabgp,'Title','Spectrum');
    tab2 =  uitab(tabgp,'Title',['MS',char(185),' Detail']);
    tab3 = uitab(tabgp,'Title',['MS',char(178),' Detail']);
    
    spectraplot = axes('parent',tab1,'Units','normalized','OuterPosition',[0 0 .85 .45],...
        'TickDir','out','XLim',[0,2000]);
    errorplot = axes('parent',tab1,'Units','normalized','OuterPosition',[0 .45 .85 .1],...
        'XLim',[0,2000]);
    textsummary = uitable('parent',tab1,'Units','normalized','OuterPosition',[.6 .6 .4 .4],...
        'RowName',[]);
    plotiontype = uicontrol('parent',tab1,'Style','checkbox','String','Show ion type',...
        'Units','normalized','OuterPosition',[.8 .52 .125 .03],...
        'Value',1,'Callback',@plotiontype_Callback);
    plotmz = uicontrol('parent',tab1,'Style','checkbox','String','m/z','Units','normalized',...
        'OuterPosition',[.8 .49 .13 .03],'Value',1,...
        'Callback',@plotmz_Callback);
    plotcharge = uicontrol('parent',tab1,'Style','checkbox','String','Charge',...
        'Units','normalized','OuterPosition',[.8 .46 .08 .03],'Value',1,...
        'Callback',@plotcharge_Callback);
    plotcomp = uicontrol('parent',tab1,'Style','checkbox','String','Composition',...
        'Units','normalized','OuterPosition',[.8 .43 .15 .03],'Value',1,...
        'Callback',@plotcomp_Callback);
    plottoppknum = uicontrol('parent',tab1,'Style','edit','Units','normalized',...
        'OuterPosition',[.85 .4 .03 .03],'String',num2str(plottpnum));
    foundpeaktable_handle = uitable('parent',tab3,'Units','normalized',...
        'OuterPosition',[.01 .01 .98 .48],'RearrangeableColumns','on');
    chkboxstat = [plotmz.Value,plotiontype.Value,plotcomp.Value,plotcharge.Value];
    plottoppknumtext1 = uicontrol('parent',tab1,'Style','text','String','Plot top',...
        'Units','normalized','OuterPosition',[.8 .395 .05 .03],...
        'HorizontalAlignment','right');
    plottoppknumtext2 = uicontrol('parent',tab1,'Style','text','String','peaks.',...
        'Units','normalized','OuterPosition',[.88 .395 .05 .03],...
        'HorizontalAlignment','left');
    retext = uicontrol('parent',tab1,'Style','pushbutton','String','Update',...
        'Units','normalized','OuterPosition',[.93 .43 .07 .07],'Value',1,...
        'Callback',@retext_Callback);
    
    thisTICplot = axes('Position',[.137 .407 .76 .216],'parent',tab2);
    monoisoplot = axes('Position',[.137 .74 .62 .216],'parent',tab2);
    totTICplot = axes('Position',[.137 .078 .76 .216],'parent',tab2);
    % Time, MonoIsotopic Plot: "tmip"
    tmipval = uicontrol('parent',tab2,'Style','text','String',['+/- ',num2str(tmipwindow),...
        ' min(s)'],'Units','normalized','Position',[.8 .79 .2 .1],...
        'FontSize',16);
    tmipinc = uicontrol('parent',tab2,'Style','pushbutton','String','+','Units','normalized',...
        'Position',[.8 .92 .07 .07],'Value',1,...
        'FontSize',24,'Callback',@tmipinc_Callback);
    tmipred = uicontrol('parent',tab2,'Style','pushbutton','String','-','Units','normalized',...
        'Position',[.8 .74 .07 .07],'Value',1,...
        'FontSize',24,'Callback',@tmipred_Callback);
    
    %% plot spectra
    peakMatchIndex_merge=~cellfun(@isempty,matchresult.peakmatchindex);
    [isoIndex,~]=findIsotopePeaks(spectrum,MS2tol,MS2tolUnit);
    for i=1:length(isoIndex)
        if any(peakMatchIndex_merge(isoIndex.mono(i)))
            peakMatchIndex_merge(isoIndex(i).iso,1)=1;
            peakMatchIndex_merge(isoIndex(i).iso,2)=1;
            peakMatchIndex_merge(isoIndex(i).iso,3)=1;
        else
            peakMatchIndex_merge(isoIndex(i).iso,1)=0;
            peakMatchIndex_merge(isoIndex(i).iso,2)=0;
            peakMatchIndex_merge(isoIndex(i).iso,3)=0;
        end
    end
    width1=0.01;
    
    % ID'd pks
    hold(spectraplot,'on');
    % unmatched: black
    bar(spectraplot,spectrum(sum(peakMatchIndex_merge,2) == 0,1),...
        spectrum(sum(peakMatchIndex_merge,2) == 0,2),...
        width1,'k','EdgeColor','k','FaceColor','k','BarWidth',0.01);
    % b/B-ion: blue
    bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,3)),1),...
        spectrum(logical(peakMatchIndex_merge(:,3)),2),...
        width1,'b','EdgeColor','b','FaceColor','b','BarWidth',0.01);
    % y/Y-ion: red
    bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,2)),1),...
        spectrum(logical(peakMatchIndex_merge(:,2)),2),...
        width1,'r','EdgeColor','r','FaceColor','r','BarWidth',0.01);
    % all others: green
    bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,1)),1),...
        spectrum(logical(peakMatchIndex_merge(:,1)),2),...
        width1,'FaceColor',[30,106,21]/128,'EdgeColor',[30,106,21]/128,...
        'BarWidth',0.01);
    
    legendplot = axes('parent',tab1,'Units','normalized','Position',[.85 0 .15 .35]);
    hold(legendplot,'on');
    axis(legendplot,'equal')
    set(legendplot,'XLim',[0,5]);
    axis(legendplot,'off')
    set(get(spectraplot,'XLabel'),'String','mz');
    set(spectraplot,'FontSize',12,'fontweight','bold');
    set(get(spectraplot,'YLabel'),'String','Ion Current');
    ylim =get(spectraplot,'ylim');
    set(spectraplot,'ylim',[ylim(1) ylim(2)*1.16]);
    set(spectraplot,'ClippingStyle','rectangle')
    foundPeaks_relativeintensity = num2cell([foundPeaks.Intensity]/max([foundPeaks.Intensity])*100);
    [foundPeaks.RelativeIntensity] = foundPeaks_relativeintensity{:};
    unichg = unique([foundPeaks.charge]);
    vertexes = [.01,.01;.99,.01;.99,.99;.01,.99];
    for i = 1:length(unichg)
        patch(vertexes(:,1),vertexes(:,2)-1+i,textcolor(unichg(i)));
        text(1.2,-.5+i,['charge ',num2str(unichg(i))]);
    end
    
    monosaclist = Glycan.gly1let;
    [unimzTheo,~,unimzTheoind] = unique([foundPeaks.mzTheo]);
    annotatetext = cell(length(unimzTheo),5);
    errorplotxy = zeros(length(unimzTheo),3);
    allpid = [foundPeaks.peakIndex];
    for i = 1 : length(unimzTheo)
        thismzind = unimzTheoind == i;
        thismzind_first = find(thismzind,1);
        corrpid = allpid(thismzind);
        corrspec = spectrum(corrpid,:);
        [locationy,locmaxind] = max(corrspec(:,2));
        locationx = double(corrspec(locmaxind,1));
        locationy = double(locationy);
        charge = foundPeaks(thismzind_first).charge;
        tempannotatetext = foundPeaks(thismzind_first).iontype;
        if strcmpi(tempannotatetext,'none') || strcmpi(tempannotatetext,'ptm')
            % INTENTIONALLY LEFT BLANK
        elseif strcmpi(tempannotatetext(1),'-')  % B/Y
            if strcmpi(tempannotatetext(2),'b')
                tempannotatetext = ['B',num2str(size(foundPeaks(thismzind_first).unitindex{2},1))];
            elseif strcmpi(tempannotatetext(2),'y')
                tempannotatetext = ['Y',num2str(size(foundPeaks(thismzind_first).unitindex{2},1))];
            end
        else  % pep break
            tempannotatetext = strsplit(tempannotatetext,'-');
            tempannotatetext = tempannotatetext{1};
            % INTENTIONALLY LEFT BLANK
        end
        
        annotatetext{i,2} = tempannotatetext;
        locationxstr = num2str(round(locationx,2));
        [~,glyMat,~] = breakGlyPep(foundPeaks(thismzind_first).sgp);
        if ~isempty(glyMat)
            thisglyms = glyMat.struct;
        else
            thisglyms = '';
        end
        thisglymslist = '';
        for j = 1:length(monosaclist)
            tempmsnum = length(strfind(thisglyms,monosaclist{j}));
            if tempmsnum > 0
                thisglymslist = [thisglymslist,monosaclist{j},'_{',num2str(tempmsnum),'}'];
            end
        end
        if ~isempty(thisglymslist)
            thisglymslist = [' (',thisglymslist,')'];
        end
        annotatetext{i,1} = locationxstr;
        annotatetext{i,3} = thisglymslist;
        annotatetext{i,4} = [num2str(charge),'^{+}'];
        annotatetext{i,5} = [locationx,locationy*1.1];
        localerrda = [foundPeaks(thismzind).DaError];
        localerrppm = [foundPeaks(thismzind).ppmError];
        [~,mindaind] = min(abs(localerrda));
        [~,minppmind] = min(abs(localerrppm));
        errorplotxy(i,:) = [locationx,localerrda(mindaind),localerrppm(minppmind)];
    end
    tempannoxy = cell2mat(annotatetext(:,5));
    [~,ind] = sort(tempannoxy(:,2),'descend');
    annotatetext = annotatetext(ind,:);
    
    
    %% 1st annotated spectrum
    % ind = min(plottpnum,size(annotatetext,1));
    
    textshown = 0;
    showhich = 1;
    hitmz = [];
    while textshown <= plottpnum && showhich <= size(annotatetext,1)
        locationx = annotatetext{showhich,5}(1);
        if ~any(hitmz + spacing > locationx & hitmz - spacing < locationx )
            charge = str2double(annotatetext{showhich,4}(1:end-4));
            text(locationx,annotatetext{showhich,5}(2),...
                [annotatetext{showhich,1},' ',annotatetext{showhich,2},' ',annotatetext{showhich,3},' ',...
                annotatetext{showhich,4}],'fontsize',8,'Clipping','on','parent',spectraplot,...
                'rotation',90,'Clipping','on','color',textcolor(charge));
            hitmz = [hitmz;locationx];
            textshown = textshown + 1;
        end
        showhich = showhich + 1;
    end
    %% error plot
    if strcmpi(MS2tolUnit,'DA')
        scatter(errorplot,errorplotxy(:,1),errorplotxy(:,2),'k');
        hold(errorplot,'on')
        set(errorplot,'FontSize',12,'fontweight','bold','YLim',[-MS2tol,MS2tol]);
        set(get(errorplot,'YLabel'),'String','Error (Da)');
    elseif strcmpi(MS2tolUnit,'PPM')
        scatter(errorplot,errorplotxy(:,1),errorplotxy(:,3),'k');
        hold(errorplot,'on')
        set(errorplot,'FontSize',12,'fontweight','bold','YLim',[-MS2tol,MS2tol]);
        set(get(errorplot,'YLabel'),'String','Error (ppm)');
    end
    
    colnames = fieldnames(foundPeaks);
    [~,newindex] = sortrows({foundPeaks.iontype2}');
    keepcolindex =[];
    for i = 1 : length(colnames)
        if strcmpi(colnames{i},'charge')  % 1
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'mass')  % 2
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'mzTheo')  % 3
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'mzExpt')  % 4
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(MS2tolUnit,'da') && strcmpi(colnames{i},'DaError')  % 5
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(MS2tolUnit,'ppm') && strcmpi(colnames{i},'ppmError')  % 5
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'iontype')  % 6
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'sgp')  % 7
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'iontype2')  % 8
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'Intensity')  % 9
            keepcolindex = [keepcolindex,i];
        elseif strcmpi(colnames{i},'RelativeIntensity')  % 10
            keepcolindex = [keepcolindex,i];
        end
    end
    colnames = colnames(keepcolindex);
    scoretabledata = cell(length(newindex),length(colnames) + 1);
    for i = 1: length(newindex)
        for j = 1:length(colnames)
            if ismember(j,[2,3,4,5])  % mass, mzTheo, mzExpt, Da/ppmerror, relative intensity
                scoretabledata{i,j} = num2str(round(foundPeaks(newindex(i)).(colnames{j}),2));
            elseif ismember(j,[9,10])  % Intensity
                scoretabledata{i,j} = num2str(round(foundPeaks(newindex(i)).(colnames{j}),0));
            else
                scoretabledata{i,j} = foundPeaks(newindex(i)).(colnames{j});
            end
        end
        nmpg = [num2str(foundPeaks(newindex(i)).nmFrag),' ',...
            num2str(foundPeaks(newindex(i)).npFrag),' ',...
            num2str(foundPeaks(newindex(i)).ngFrag)];  %nmFrag/npFrag/ngFrag
        scoretabledata{i,j+1} = nmpg;
    end
    colnames = [colnames;'nm/np/ngFrag'];
    % charge  mass  mzTheo  mzExpt  Error  iontype
    % sgp  iontype2  Intensity  RelaIntensity  nm/np/ngFrag
    columnwidth = {40,55,55,55,55,200,...
        240,80,55,100,100};
    columnformat = {'numeric','longg','longg','longg','longg','char',...
        'char','char','numeric','numeric','char'};
    set(foundpeaktable_handle,'data',scoretabledata,...
        'ColumnFormat',columnformat);
    set(foundpeaktable_handle,'ColumnWidth',columnwidth);
    set(foundpeaktable_handle,'ColumnName',colnames);
end
%% 2 x TICs
hold(thisTICplot,'on');
hold(totTICplot,'on');
totTICdata = [allretime,cell2mat(alltotIonCurrent)];
totTICdata = totTICdata(allmslvl == 1,:);  % chroma plot
bar(thisTICplot,thisTICdata(:,1),thisTICdata(:,2),'Edgecolor','k','FaceColor','k');
bar(totTICplot,totTICdata(:,1),totTICdata(:,2),'Edgecolor','k','FaceColor','k');
thisticylim = thisTICplot.YLim;
plot(thisTICplot,[parentretime,parentretime],[0,max(thisticylim)],'r');
totticylim = totTICplot.YLim;
plot(totTICplot,[parentretime,parentretime],[0,max(totticylim)],'r');
thisTICplot.XLabel.String = 'Time (min)';
thisTICplot.YLabel.String = 'Precursor Ion Current';
totTICplot.XLabel.String = 'Time (min)';
totTICplot.YLabel.String = 'Total Ion Current';
% xlim = get(totTICplot,'XLim');
% set(thisTICplot,'XLim',xlim);

%% AUC
[~,thisretimeind] = min(abs(thisTICdata(:,1) - parentretime));
leftrtind = thisretimeind;
rightrtind = thisretimeind;
go_on = true;
while go_on
    if leftrtind >= 1 && thisTICdata(leftrtind,2) > 0
        leftrtind = leftrtind - 1;
    else
        go_on = false;
    end
end
go_on = true;
while go_on
    if rightrtind <= size(thisTICdata,1) && thisTICdata(rightrtind,2) > 0
        rightrtind = rightrtind + 1;
    else
        go_on = false;
    end
end
leftrtind = max(leftrtind,1);
rightrtind = min(rightrtind,size(thisTICdata,1));
if leftrtind ~= rightrtind
    thisAUC = trapz(thisTICdata(leftrtind:rightrtind,1),thisTICdata(leftrtind:rightrtind,2));
else
    thisAUC = thisTICdata(leftrtind,2);
end
ticstrpos = get(thisTICplot,'YLim') * 0.5;
text(thisTICplot,parentretime,ticstrpos(2),num2str(round(thisAUC,0)));

plotmonoisodist(msdata,score,monoisoplot,tmipval,tmipwindow);

summarydata = {'Scan number',score.Scan;...
    'Charge',score.Charge;...
    'Fragmentation mode',score.Fragmode;...
    'SGP',score.SGP;...
    'Enscore',score.Enscore;...
    'Decoy enscore',score.DecoyEnscore;...
    '% Gly. bond fragged',round(100*score.Glycov);...
    '% Pep. bond frag(w/o gly. residue)',round(100*score.Fracpepfragmatch);...
    '% Pep. bond frag(w/ gly. residue)',round(100*score.Fracglyfragmatch);...
    'Peak AUC',thisAUC;...
    'Retention time (min.)',round(ms2retime,2);...
    'Parent scan re. time (min.)',round(parentretime,2);...
    'Parent ion mass',num2str(round(originalprecursorMz(1),4));...
    'Ion monoiso. mass (esti.)',num2str(round(precursorMz,4));...
    'Candi. theo. mass',num2str(round(score.Theo,4));...
    'Peak lag',score.PeakLag;...
    '% Theo. frag. ion match',score.PercentIonMatch;...
    'P-value',score.Pvalue;...
    'Max. top 10 peaks hit',score.Top10;...
    'NpFrag',score.NpFrag;...
    'NgFrag',score.NgFrag;...
    'NmFrag',score.NmFrag};
set(textsummary,'Data',summarydata,'ColumnName',{'Name','Value'},'ColumnWidth',{135,165});

spectraplotxlim = get(spectraplot,'XLim');
plot(errorplot,[min(spectraplotxlim),max(spectraplotxlim)],[0,0],'k');
set(errorplot,'XLim',spectraplotxlim,'TickDir','out','XTick',[],'XTickLabel',[]);

%% 2 structure figures
specialoptions = getdgfrag({p,g,m},theofrag,matchresult);
specialoptions_simp = specialoptions;
for i = 1:length(specialoptions_simp)
    tempstruct = specialoptions_simp{i};
    fldnms = {'R','NR','pepN','pepC'};
    for j = 1:length(fldnms)
        if isfield(tempstruct,fldnms{j})
            tempfld = tempstruct.(fldnms{j});
            tempfld(:,1) = cell(size(tempfld,1),1);
            tempstruct.(fldnms{j}) = tempfld;
        end
    end
    specialoptions_simp{i} = tempstruct;
end
sgpplotB = axes('Position',[.01 .5 .98 .5],'parent',tab3);  % detailed DG figure @ tab 3
hold(sgpplotB,'on');
figinfoB = drawglycan(sgp,'figurehandle',sgpplotB,'inputformat','SGP1',...
    'perpendicularmonosac',{'Deoxyhexose','Pentose','Fuc','Xyl'}...
    ,'specialoptions',specialoptions);
tcpB = figinfoB.tcp + [-.5 -.5;.5 .5];  % add margin to drawglycan figure
set(sgpplotB,'XLim',[tcpB(1,1),tcpB(2,1)],'YLim',[tcpB(1,2),tcpB(2,2)]);
axis(sgpplotB,'off')
sgpplotB.Units = 'normalized';

sgpplotA = axes('Position',[.0 .55 .6 .45],'parent',tab1);  % concise DG figure @ tab 1
hold(sgpplotA,'on');
figinfoA = drawglycan(sgp,'figurehandle',sgpplotA,'inputformat','SGP1',...
    'perpendicularmonosac',{'Deoxyhexose','Pentose','Fuc','Xyl'}...
    ,'specialoptions',specialoptions_simp);
tcpA = figinfoA.tcp + [-.5 -.5;.5 .5];  % add margin to drawglycan figure
set(sgpplotA,'XLim',[tcpA(1,1),tcpA(2,1)],'YLim',[tcpA(1,2),tcpA(2,2)]);
axis(sgpplotA,'off')
sgpplotA.Units = 'normalized';

if ~isempty(scoretablefilename)
    FID = fopen(scoretablefilename,'w');
    foundPeaks_write = rmfield(foundPeaks,'unitindex');
    struct2csv(foundPeaks_write,FID);
    fclose(FID);
end

%% Finally
if(figdisplay)
    set(h,'Visible','on');
end

    function plotiontype_Callback(source,eventdata)
        chkboxstat = [plotmz.Value,plotiontype.Value,plotcomp.Value,plotcharge.Value];
    end

    function plotmz_Callback(source,eventdata)
        chkboxstat = [plotmz.Value,plotiontype.Value,plotcomp.Value,plotcharge.Value];
    end

    function plotcharge_Callback(source,eventdata)
        chkboxstat = [plotmz.Value,plotiontype.Value,plotcomp.Value,plotcharge.Value];
    end

    function plotcomp_Callback(source,eventdata)
        chkboxstat = [plotmz.Value,plotiontype.Value,plotcomp.Value,plotcharge.Value];
    end

    function retext_Callback(source,eventdata)
        cla(spectraplot);
        hold(spectraplot,'on');
        bar(spectraplot,spectrum(sum(peakMatchIndex_merge,2) == 0,1),...
            spectrum(sum(peakMatchIndex_merge,2) == 0,2),width1,'FaceColor','k','EdgeColor','k');
        bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,3)),1),...
            spectrum(logical(peakMatchIndex_merge(:,3)),2),width1,'FaceColor','b','EdgeColor','b');
        bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,2)),1),...
            spectrum(logical(peakMatchIndex_merge(:,2)),2),width1,'FaceColor','k','EdgeColor','r');
        bar(spectraplot,spectrum(logical(peakMatchIndex_merge(:,1)),1),...
            spectrum(logical(peakMatchIndex_merge(:,1)),2),...
            width1,'FaceColor',[30,106,21]/128,'EdgeColor',[30,106,21]/128);
        plottpnum = str2double(plottoppknum.String);
        textshown = 0;
        showhich = 1;
        hitmz = [];
        while textshown <= plottpnum && showhich <= size(annotatetext,1)
            tlocationx = annotatetext{showhich,5}(1);
            if ~any(hitmz + spacing > tlocationx & hitmz - spacing < tlocationx )
                tatext = '';
                if chkboxstat(1)
                    tatext = [tatext,' ',annotatetext{showhich,1}];
                end
                if chkboxstat(2)
                    tatext = [tatext,' ',annotatetext{showhich,2}];
                end
                if chkboxstat(3)
                    tatext = [tatext,' ',annotatetext{showhich,3}];
                end
                if chkboxstat(4)
                    tatext = [tatext,' ',annotatetext{showhich,4}];
                end
                tcharge = str2double(annotatetext{showhich,4}(1:end-4));
                text(tlocationx,annotatetext{showhich,5}(2),tatext,'fontsize',8,'Clipping','on',...
                    'parent',spectraplot,'rotation',90,'Clipping','on','color',textcolor(tcharge));
                hitmz = [hitmz;locationx];
                textshown = textshown + 1;
            end
            showhich = showhich + 1;
        end
    end

    function tmipinc_Callback(source,eventdata)
        if tmipwindow == 0.25 || tmipwindow == 0.5
            tmipwindow = tmipwindow * 2;
        elseif tmipwindow == 0
            tmipwindow = 0.25;
        elseif tmipwindow == 10
            % INTENTIONALLY LET BLANK
        else
            tmipwindow = tmipwindow + 1;
        end
        plotmonoisodist(msdata,score,monoisoplot,tmipval,tmipwindow);
    end

    function tmipred_Callback(source,eventdata)
        if tmipwindow == 1 || tmipwindow == 0.5
            tmipwindow = tmipwindow / 2;
        elseif tmipwindow == 0.25
            tmipwindow = 0;
        elseif tmipwindow == 0
            tmipwindow = 0;
        else
            tmipwindow = tmipwindow - 1;
        end
        plotmonoisodist(msdata,score,monoisoplot,tmipval,tmipwindow);
    end
end

function foundPeaks = easymatch(spectrum,sgp,result,theofrag)
% Fix compatibility issues with old code, also to get theo-expt error
% values.
foundPeaks = [];
pkmatchind = result.peakmatchindex;
fcount = 1;
totalmonosac = 0;
intaction = theofrag(ismember({theofrag.type},'none'));
totalmonosac = length(intaction.unitindex{2});
for i = 1:length(pkmatchind)
    temppmi = pkmatchind{i};
    if ~isempty(temppmi)
        for j = 1:size(temppmi,1)
            foundPeaks(fcount).charge=temppmi(j,2);
            foundPeaks(fcount).mass=theofrag(temppmi(j,1)).mz - 1.007825032;
            foundPeaks(fcount).mzExpt=spectrum(i,1);
            foundPeaks(fcount).mzTheo=(theofrag(temppmi(j,1)).mz - 1.007825032)/temppmi(j,2)...
                + 1.007825032;
            foundPeaks(fcount).DaError = spectrum(i,1)-foundPeaks(fcount).mzTheo;
            foundPeaks(fcount).ppmError=(spectrum(i,1)-foundPeaks(fcount).mzTheo)/...
                foundPeaks(fcount).mzTheo*1e6;
            iontype = theofrag(temppmi(j,1)).type;
            foundPeaks(fcount).iontype = iontype;
            foundPeaks(fcount).sgp=theofrag(temppmi(j,1)).sgp;
            if strcmpi(iontype(1),'-')
                foundPeaks(fcount).iontype2 = [upper(iontype(2)),...
                    num2str(length(theofrag(temppmi(j,1)).unitindex{2}))];
            elseif strcmpi(iontype,'none')
                foundPeaks(fcount).iontype2 = 'NONE';
            else
                nummonosac = length(theofrag(temppmi(j,1)).unitindex{2});
                if nummonosac == totalmonosac
                    foundPeaks(fcount).iontype2 = iontype;
                else
                    if any(strfind(iontype,'-'))
                        iontype = strsplit(iontype,'-');
                        foundPeaks(fcount).iontype2 = [iontype{1},' - Y',...
                            num2str(nummonosac)];
                    else
                        foundPeaks(fcount).iontype2 = iontype;
                    end
                end
            end
            foundPeaks(fcount).Intensity=spectrum(i,2);
            foundPeaks(fcount).original=sgp;
            foundPeaks(fcount).nmFrag=theofrag(temppmi(j,1)).nmFrag;
            foundPeaks(fcount).npFrag=theofrag(temppmi(j,1)).npFrag;
            foundPeaks(fcount).ngFrag=theofrag(temppmi(j,1)).ngFrag;
            foundPeaks(fcount).peakIndex = i;
            foundPeaks(fcount).unitindex = theofrag(temppmi(j,1)).unitindex;
            fcount = fcount + 1;
        end
    end
end
end