function [thisfigsto,theofragsto,annobuttonsto,CombiES,DecoyCombiES,...
    structfigsto,currentind,spectrum,theores] = ...
    createtriggeredannotations(msdata,scoredata,...
    current_displaydataind,currentisomerind,scoreintdata,tgA,uiscoretable,...
    theofragsto,thisfigsto,annobuttonsto,structfigsto)
ms2tol = scoreintdata.scoreoptions.ms2tol;
ms2tolunit = scoreintdata.scoreoptions.ms2tolunit;
fragmethods = scoreintdata.scoreoptions.analyzefragmode;
if scoreintdata.scoreoptions.isHCDtrigger
    [~,ind1] = sort(scoreintdata.colnames);
    [~,ind2] = sort(fragmethods);
    [~,ind3] = sort(ind1);
    ind = ind2(ind3);
    ms2tol = ms2tol(ind);
    ms2tolunit = ms2tolunit(ind);
    fragmethods = fragmethods(ind);
end
currentind = zeros(size(current_displaydataind));
for ii = 1:length(current_displaydataind)
    currentind(ii) = current_displaydataind{ii}(currentisomerind);
end
[CombiES,DecoyCombiES] = calculateCombiES(scoredata,currentind,...
    fragmethods);
spectra = msdata.spectra;
scannum = msdata.scannum;
charge = msdata.charge;
current_scoredata= scoredata(currentind(currentind > 0));
numannofig = length(currentind);
allresultscannum = [current_scoredata.Scan];
if contains(scoreintdata.sliminput.pepfile,'_O_glycopep_digest.txt')
    uiscoretable_colnames = {'Scan';'Fragmode';'Theo';'Enscore';'Glycov';...
        'Fracpepfragmatch';'Fracglyfragmatch';'PercentIonMatch';'Top10';'Pvalue';...
        'PeakLag';'HtCenter';'HtAvg';'Y0Y1Y2';'Mono';...
        'Expt';'DecoyRatio';'SelectPeak';'NpFrag';'NgFrag';...
        'NmFrag';'DecoyEnscore';'Retime';'ProteinID';'MDiff';'localization'};
else
        uiscoretable_colnames = {'Scan';'Fragmode';'Theo';'Enscore';'Glycov';...
        'Fracpepfragmatch';'Fracglyfragmatch';'PercentIonMatch';'Top10';'Pvalue';...
        'PeakLag';'HtCenter';'HtAvg';'Y0Y1Y2';'Mono';...
        'Expt';'DecoyRatio';'SelectPeak';'NpFrag';'NgFrag';...
        'NmFrag';'DecoyEnscore';'Retime';'ProteinID';'MDiff'};
end
uiscoretable_dispdata = cell(length(current_scoredata),length(uiscoretable_colnames));
for ii = 1:length(uiscoretable_colnames)
    uiscoretable_dispdata(:,ii) = {current_scoredata.(uiscoretable_colnames{ii})};
end
set(uiscoretable,'ColumnName',uiscoretable_colnames,'Data',uiscoretable_dispdata);
figposfactor = 0.98/numannofig;
sgp = current_scoredata(1).SGP;
theores = cell(numannofig,1);
theofrag = cell(numannofig,1);
spectrum = cell(numannofig,1);

if isempty(thisfigsto)
    thisfigsto = cell(1,numannofig);
    annobuttonsto = cell(1,numannofig);
    structfigsto = cell(1,numannofig);
end

if isempty(theofragsto)
    isomerind = zeros(length(current_displaydataind{1}),1);
    for ii = 1:length(isomerind)
        for jj = 1:length(current_displaydataind)
            tempind = current_displaydataind{jj}(ii);
            if tempind > 0
                isomerind(ii) = tempind;
            end
        end
    end
    isomersgpcomp = {scoredata(isomerind).ProteinID};
    isomersgpseq = {scoredata(isomerind).SGP};
    scoreintdata.isomersgpseq = isomersgpseq;
    theofragsto = createtheofragsto(isomersgpcomp,fragmethods,scoreintdata);
end

denoisingoptions.fragmode = scoreintdata.scoreoptions.analyzefragmode;
denoisingoptions.cutoffmed = scoreintdata.scoreoptions.cutoffmed;
denoisingoptions.fracmax = scoreintdata.scoreoptions.fracmax;
for ii = 1:numannofig
    set(annobuttonsto{ii},'String','Annotate');
    if currentind(ii) == 0 && ~isempty(annobuttonsto{ii})
        set(annobuttonsto{ii},'visible','off');
        set(thisfigsto{ii},'visible','off');
        cla(thisfigsto{ii});
        cla(structfigsto{ii});
    else
        if currentind(ii) > 0
            selectedres = scoredata(currentind(ii));
            if isempty(thisfigsto{ii})
                thisannofig = axes('Parent',tgA,'Unit','normalized',...
                    'Position',[0.05+figposfactor*(ii-1),.07,figposfactor*.85,.35]);
                title(thisannofig,selectedres(1).Fragmode);
                xlabel(thisannofig,'m/z');
                ylabel(thisannofig,'Relative Intenisty')
                hold(thisannofig,'on');
                thisfigsto{ii} = thisannofig;
                anno = uicontrol('Parent',tgA,'Style','pushbutton','String',...
                    'Annotate','Unit','normalized','visible','on',...
                    'Position',[.05+figposfactor*(ii-.7),.46,figposfactor/3,.03],...
                    'Callback',{@annobutton_callback,ii});
                annobuttonsto{ii} = anno;
                thisstructfig = axes('Parent',tgA,'Unit','normalized',...
                    'Position',[0.32+figposfactor*0.68*(ii-1),.5,figposfactor*.51,.3]);
                axis(thisstructfig,'off')
                structfigsto{ii} = thisstructfig;
            end
            thisscannum = selectedres.Scan;
            thisannofig = thisfigsto{ii};
            thisstructfig = structfigsto{ii};
            
            cla(thisannofig);
            set(annobuttonsto{ii},'visible','on');
            cla(thisstructfig);
            
            tempspec = spectra{scannum == thisscannum};
            thischarge = charge(scannum == thisscannum);
            fragmethod = upper(selectedres(1).Fragmode);
            thisspec = spectradenoising(tempspec,fragmethod,thischarge,...
                glypepMW(sgp),ms2tol(ii),ms2tolunit{ii},denoisingoptions);
            spectrum{ii} = thisspec;
            theofrag{ii} = theofragsto{currentisomerind,ii};
            cisoptions.maxlag = scoreintdata.scoreoptions.maxlag;
            cisoptions.selectpeak = scoreintdata.scoreoptions.selectpeak;
            theores{ii} = calcithscore(thisspec,[theofrag{ii}.mz],thischarge,...
                ms2tol(ii),ms2tolunit{ii},2,cisoptions);
            matchind = ~cellfun(@isempty,theores{ii}.peakmatchindex);
            bar(thisannofig,thisspec(~matchind,1),thisspec(~matchind,2)/max(thisspec(:,2))*100,...
                'edgecolor','r','FaceColor','r','BarWidth',0.01);
            bar(thisannofig,thisspec(matchind,1),thisspec(matchind,2)/max(thisspec(:,2))*100,...
                'edgecolor','g','FaceColor','g','BarWidth',0.01);
            for jj = 1:size(theores{ii}.peakmatchindex,1)
                for kk = 1:size(theores{ii}.peakmatchindex{jj},1)
                    tempmatchedtheofrag = theofrag{ii}(theores{ii}.peakmatchindex{jj}(kk,1));
                    if ~isequal([tempmatchedtheofrag.nmFrag,tempmatchedtheofrag.ngFrag],[0,0]) && tempmatchedtheofrag.npFrag > 0
                        text(thisannofig,thisspec(jj,1),thisspec(jj,2)/max(thisspec(:,2))*100*1.1,'*');
                    end
                end
            end
            if isfield(theofrag{ii},'uniqueness')
                pkmatchind = theores{ii}.peakmatchindex;
                for jj = 1:length(pkmatchind)
                    temppkmi = pkmatchind{jj};
                    if ~isempty(temppkmi)
                        uniquenessstr = '';
                        for kk = 1:size(temppkmi,1)
                            frag = theofrag{ii}(temppkmi(kk,1));
                            fraguniqueness = frag.uniqueness;
                            if any(strfind(fraguniqueness,'unique')) || any(strfind(fraguniqueness,'shared'))
                                uniquenessstr = [uniquenessstr,fraguniqueness,', '];
                            end
                        end
                        if ~isempty(uniquenessstr)
                            txtpos = thisspec(jj,:);
                            text(thisannofig,txtpos(1),txtpos(2)/max(thisspec(:,2))*100*1.1,uniquenessstr,'rotation',90);
                        end
                    end
                end
            end
            set(thisannofig,'visible','on');
            thisfigsto{ii} = thisannofig;
            set(annobuttonsto{ii},'Callback',{@annobutton_callback,ii});
            
            [p,g,m] = breakGlyPep(sgp);
            specialoptions = getdgfrag({p,g,m},theofrag{ii},theores{ii});
            specialoptions_simp = specialoptions;
            for jj = 1:length(specialoptions_simp)
                tempstruct = specialoptions_simp{jj};
                fldnms = {'R','NR','pepN','pepC'};
                for kk = 1:length(fldnms)
                    if isfield(tempstruct,fldnms{kk})
                        tempfld = tempstruct.(fldnms{kk});
                        if strcmpi(fldnms{kk},'pepN') || strcmpi(fldnms{kk},'pepC')
                            for mm = 1:size(tempfld,1)
                                if any(strfind(tempfld{mm,1},'*'))
                                    tempfld{mm,1} = '*';
                                end
                                firstnumpos = -1;
                                for nn = length(tempfld{mm,1}):-1:1
                                    if double(tempfld{mm,1}(nn)) <= 57 &&  double(tempfld{mm,1}(nn)) >= 48
                                        firstnumpos = nn;
                                    end
                                end
                                if firstnumpos > 0
                                    tempfld{mm,1} = tempfld{mm,1}(1:firstnumpos - 1);
                                end
                            end
                        elseif strcmpi(fldnms{kk},'R') || strcmpi(fldnms{kk},'NR')
                            tempfld(:,1) = cell(size(tempfld,1),1);
                        end
                        tempstruct.(fldnms{kk}) = tempfld;
                    end
                end
                specialoptions_simp{jj} = tempstruct;
            end
            hold(thisstructfig,'on');
            figinfoA = drawglycan(sgp,'figurehandle',thisstructfig,'inputformat','SGP1',...
                'perpendicularmonosac',{'Deoxyhexose','Pentose','Fuc','Xyl'}...
                ,'specialoptions',specialoptions_simp,'fontsize',8);
            tcpA = figinfoA.tcp + [-.5 -.5;.5 .5];  % add margin to drawglycan figure
            set(thisstructfig,'XLim',[tcpA(1,1),tcpA(2,1)],'YLim',[tcpA(1,2),tcpA(2,2)]);
            axis(thisstructfig,'off')
            thisstructfig.Units = 'normalized';
        end
    end
end
    function annobutton_callback(src,event,x)
        genannofig(allresultscannum(x),spectrum{x},current_scoredata(x),...
            ms2tol(x),ms2tolunit{x},...
            msdata,theores{x},theofrag{x},'');
    end
end