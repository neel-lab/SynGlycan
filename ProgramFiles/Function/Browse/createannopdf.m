function createannopdf(msdata,scoredata,displaydataind,displayoptions)
% CREATEANNOPDF: create annotationg
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
% Date Lastly Updated: 11/10/2020

answer = questdlg('This operation will export all data displayed in table. Continue?', ...
    'Confirm Result Data Range', ...
    'Yes, proceed','No','Yes, proceed');
switch answer
    case 'Yes, proceed'
        spectra = msdata.spectra;
        scannum = msdata.scannum;
        charge = msdata.charge;
        savetofolder = uigetdir('','Select Output Folder');
        scoreintdata = displayoptions.scoreintdata;
        denoisingoptions.fragmode = scoreintdata.scoreoptions.analyzefragmode;
        denoisingoptions.cutoffmed = scoreintdata.scoreoptions.cutoffmed;
        denoisingoptions.fracmax = scoreintdata.scoreoptions.fracmax;
        if ischar(savetofolder)
            if iscell(displaydataind)
                outputdisplaydataind = displaydataind(:);
                outputdisplaydataind = outputdisplaydataind(~cellfun(@isempty,outputdisplaydataind));
                tempoutputdisplaydataind = [];
                for ii = 1:length(outputdisplaydataind)
                    temp2outputdisplaydataind = outputdisplaydataind{ii};
                    tempoutputdisplaydataind = [tempoutputdisplaydataind;temp2outputdisplaydataind(:)];
                end
                outputdisplaydataind = tempoutputdisplaydataind;
            else
                outputdisplaydataind = displaydataind(:);
                outputdisplaydataind(outputdisplaydataind == 0) = [];
            end
            result = scoredata(outputdisplaydataind);
            textcolor = 'kbrmgyc';
            for ii = 1:length(result)
                sgp = result(ii).SGP;
                [p,g,m] = breakGlyPep(sgp);
                fragmethod = upper(result(ii).Fragmode);
                scannumind = scannum == result(ii).Scan;
                fragmethodind = ismember(upper(scoreintdata.sliminput.fragmode),fragmethod);
                ms2tol = scoreintdata.sliminput.ms2tol(fragmethodind);
                ms2tolunit = scoreintdata.sliminput.ms2tolunit{fragmethodind};
                theofrag = createtheofragsto({result(ii).ProteinID},{result(ii).Fragmode},scoreintdata);
                theofrag = theofrag{1};
                fig = figure('Units','pixels','Position',[0,0,800,360],'color','none','Visible',false);
                set(fig,'Units','normalized');
                header = uicontrol('Parent',fig,'style','text','Units','normalized',...
                    'Position',[.02,.94,.65,.04],'backgroundcolor','white');
                subheader = uicontrol('Parent',fig,'style','text','Units','normalized',...
                    'Position',[.02,.91,.7,.07],'backgroundcolor','white');
                anno{3} = axes('Parent',fig,'Units','normalized',...
                    'OuterPosition',[.02,.02,.7,.9]);
                frag{3} = axes('Parent',fig,'Units','normalized',...
                    'Position',[.67,.02,.31,.9],'Visible',false);
                set(header,'String',sgp);
                tempspec = spectra{scannumind};
                thischarge = charge(scannumind);
                thisspec = spectradenoising(tempspec,fragmethod,thischarge,...
                    glypepMW(sgp),ms2tol,ms2tolunit,denoisingoptions);
                cisoptions.maxlag = scoreintdata.sliminput.maxlag;
                cisoptions.selectpeak = [];
                tempres = calcithscore(thisspec,[theofrag.mz],thischarge,...
                    ms2tol,ms2tolunit,2,cisoptions);
                peakMatchIndex_merge=~cellfun(@isempty,tempres.peakmatchindex);
                [isoIndex,~]=findIsotopePeaks(thisspec,ms2tol,ms2tolunit);
                for jj=1:length(isoIndex)
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
                % ID'd pks
                hold(anno{3},'on');
                % unmatched: black
                barspec = thisspec;
                barspec(sum(peakMatchIndex_merge,2) ~= 0,2) = 0;
                bar(anno{3},barspec(:,1),barspec(:,2),...
                    width1,'k','EdgeColor','k','FaceColor','k');
                % b/B-ion: blue
                barspec = thisspec;
                barspec(~logical(peakMatchIndex_merge(:,3)),2) = 0;
                bar(anno{3},barspec(:,1),barspec(:,2),...
                    width1,'b','EdgeColor','b','FaceColor','b');
                % y/Y-ion: red
                barspec = thisspec;
                barspec(~logical(peakMatchIndex_merge(:,2)),2) = 0;
                bar(anno{3},barspec(:,1),barspec(:,2),...
                    width1,'r','EdgeColor','r','FaceColor','r');
                % all others: green
                barspec = thisspec;
                barspec(~logical(peakMatchIndex_merge(:,1)),2) = 0;
                bar(anno{3},barspec(:,1),barspec(:,2),...
                    width1,'FaceColor',[30,106,21]/128,'EdgeColor',[30,106,21]/128);
                set(get(anno{3},'XLabel'),'String','mz');
                set(anno{3},'FontSize',12,'fontweight','bold');
                set(get(anno{3},'YLabel'),'String','Ion Current');
                ylim =get(anno{3},'ylim');
                set(anno{3},'ylim',[ylim(1) ylim(2)*1.16]);
                set(anno{3},'ClippingStyle','rectangle');
                foundPeaks = easymatch(thisspec,sgp,tempres,theofrag);
                [unimzTheo,~,unimzTheoind] = unique([foundPeaks.mzTheo]);
                annotatetext = cell(length(unimzTheo),5);
                allpid = [foundPeaks.peakIndex];
                for jj = 1:length(unimzTheo)
                    thismzind = unimzTheoind == jj;
                    thismzind_first = find(thismzind,1);
                    corrpid = allpid(thismzind);
                    corrspec = thisspec(corrpid,:);
                    [locationy,locmaxind] = max(corrspec(:,2));
                    locationx = double(corrspec(locmaxind,1));
                    locationy = double(locationy);
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
                    annotatetext{jj,1} = [locationx,locationy*1.1];
                    annotatetext{jj,2} = tempannotatetext;
                    annotatetext{jj,3} = [num2str(foundPeaks(thismzind_first).charge),'^{+}'];
                end
                tempannoxy = cell2mat(annotatetext(:,3));
                [~,ind] = sort(tempannoxy(:,1),'ascend');
                annotatetext = annotatetext(ind,:);
                for jj = 1:size(annotatetext,1)
                    locationx = annotatetext{jj,1}(1);
                    pkcharge = str2double(annotatetext{jj,3}(1:end-4));
                    text(locationx,annotatetext{jj,1}(2),...
                        [annotatetext{jj,2},' ',annotatetext{jj,3}],...
                        'fontsize',8,'Clipping','on','parent',anno{3},...
                        'rotation',90,'Clipping','on','color',textcolor(pkcharge));
                end
                set(subheader,'String',['Scan = ',num2str(result(ii).Scan),...
                    ' Charge = ',num2str(thischarge),...
                    ' FragMode = ',result(ii).Fragmode,...
                    ' Enscore = ',num2str(result(ii).Enscore)],'HorizontalAlignment','Left');
                specialoptions = getdgfrag({p,g,m},theofrag,tempres);
                hold(frag{3},'on');
                drawglycan(sgp,'figurehandle',frag{3},'inputformat','SGP1',...
                    'perpendicularmonosac',{'Deoxyhexose','Pentose','Xyl','Fuc'},...
                    'fontsize',6,'aaspacing',0.4,'specialoptions',specialoptions);
                frag{3}.Units = 'normalized';
                axis(frag{3},'equal')
                set(fig,'PaperUnits','normalized',...
                    'PaperPosition',[0.02,0.02,0.96,0.96],...
                    'PaperType','usletter',...
                    'InvertHardcopy', 'off');
                scoretablefilename = fullfile(savetofolder,[num2str(ii),'.pdf']);
                print(fig,scoretablefilename,'-dpdf','-painters');
                close(fig);
            end
        end
end
end