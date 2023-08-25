function annotatespectra(msdata,scoredata,displaydataind,displayoptions,usercustom,selectedrow)

%% Set up control elements
tgA = figure('Units','Pixels','Position',[50 50 1600 728],'NumberTitle','off',...
    'CloseRequestFcn',@closereq,'KeyPressFcn',@keyboardaction);  % triggeredAnnotate / tgA
% uiscoretable_colnames = {'Scan';'Fragmode';'Theo';'Enscore';'Glycov';...
%     'Fracpepfragmatch';'Fracglyfragmatch';'PercentIonMatch';'Top10';'Pvalue';...
%     'PeakLag';'HtCenter';'HtAvg';'Y0Y1Y2';'Mono';...
%     'Expt';'DecoyRatio';'SelectPeak';'NpFrag';'NgFrag';...
%     'NmFrag';'DecoyEnscore';'Retime';'ProteinID';'MDiff'};
columnwidth = {50 50 60 40 40  ...
    70 70 70 40 50  ...
    50 50 50 150 50  ...
    50 50 50 35 35  ...
    35 50 50 50 50};
uiscoretable = uitable('Parent',tgA,'Unit','normalized','Position',[.02 .85 .78 .13],...
    'ColumnWidth',columnwidth);
prevrow = uicontrol('Parent',tgA,'Style','pushbutton','String','-Row','Unit','normalized',...
    'Position',[.89 0.95 .03 .04],'Callback',@prevrow_callback);
currentrow = uicontrol('Parent',tgA,'Style','text','String','','Unit','normalized',...
    'Position',[.92 0.945 .04 .04]);
nextrow = uicontrol('Parent',tgA,'Style','pushbutton','String','+Row','Unit','normalized',...
    'Position',[.96 0.95 .03 .04],'Callback',@nextrow_callback);
previso = uicontrol('Parent',tgA,'Style','pushbutton','String','-Iso','Unit','normalized',...
    'Position',[.89 0.9 .03 .04],'Callback',@previso_callback);
currentiso = uicontrol('Parent',tgA,'Style','text','String','','Unit','normalized',...
    'Position',[.92 0.895 .04 .04]);
nextiso = uicontrol('Parent',tgA,'Style','pushbutton','String','+Iso','Unit','normalized',...
    'Position',[.96 0.9 .03 .04],'Callback',@nextiso_callback);
jumptohelptext = uicontrol('Parent',tgA,'Style','Text','String','Jump to:','Unit','normalized',...
    'Position',[.86 0.855 .028 .03]);
jumptowhat = uicontrol('Parent',tgA,'Style','Popupmenu','String',{'Row','Scan','Isomer'},'Unit','normalized',...
    'Position',[.89 0.855 .035 .03]);
jumptowhere = uicontrol('Parent',tgA,'Style','Edit','String','','Unit','normalized',...
    'Position',[.93 0.855 .025 .03]);
rejectallisomer = uicontrol('Parent',tgA,'Style','Pushbutton','String','Reject.','Unit','normalized',...
    'Position',[.81 0.85 .03 .04],'Callback',@rejectall_callback);
jump = uicontrol('Parent',tgA,'Style','pushbutton','String','Jump!','Unit','normalized',...
    'Position',[.96 .85 .03 .04],'Callback',@jump_callback);
%ssa = uicontrol('Parent',tgA,'Style','pushbutton','String','RunSSA','Unit','normalized',...
%    'Position',[.85 .8 .032 .04],'Callback',@ssa_callback);
%suggcandi = uicontrol('Parent',tgA,'Style','pushbutton','String','MoreSGP','Unit','normalized',...
%    'Position',[.884 .8 .04 .04],'Callback',@suggcandi_callback);
%changecandidate = uicontrol('Parent',tgA,'Style','pushbutton','String','UpdateSGP','Unit','normalized',...
%    'Position',[.926 .8 .04 .04],'Callback',@changecandidate_callback);
saveres = uicontrol('Parent',tgA,'Style','pushbutton','String','Save','Unit','normalized',...
    'Position',[.968 .8 .022 .04],'Callback',@saveres_callback);
cb_specqual = uicontrol('Parent',tgA,'Style','checkbox','String','Spectrum quality','Unit','normalized',...
    'Position',[.81 0.965 .07 .02],'Value',false,'Callback',@cb_specqual_callback);
cb_ms1qual = uicontrol('Parent',tgA,'Style','checkbox','String','MS1 quality','Unit','normalized',...
    'Position',[.81 0.93 .07 .02],'Value',false,'Callback',@cb_ms1qual_callback);
cb_subjective = uicontrol('Parent',tgA,'Style','checkbox','String','User subjective','Unit','normalized',...
    'Position',[.81 0.895 .07 .02],'Value',false,'Callback',@cb_subjective_callback);
% annotategroup = uicontrol('Parent',tgA,'Style','pushbutton','String','DrawAll','Unit','normalized',...
%     'Position',[.845 .8 .037 .04],'Callback',@annotategroup_callback);
SGPwashelptext = uicontrol('Parent',tgA,'Style','text','Unit','normalized','Position',[.02 .8 .1 .03],...
    'String','Candidate SGP unchanged: ');
editbox_SGPwas = uicontrol('Parent',tgA,'Style','edit','Unit','normalized','Position',[.12 .8 .28 .03],...
    'String','');
%txt_simultaneousfrag = uicontrol('Parent',tgA,'Style','Text','String',...
%    'Default simultaneous frag. setting','Unit','normalized',...
%    'Position',[.65 .8 .15 .03]);
%togglesimufrag = uicontrol('Parent',tgA,'Style','pushbutton','String','Simu. Frag','Unit','normalized',...
%    'Position',[.81 .8 .038 .04],'Callback',@simufrag_callback);
simufragstat = false;
MS1dist = axes('Parent',tgA,'Unit','normalized','Position',[.02 .6 .25 .15]);
tmipwindow = 0.25; % 0.5 min time window for monoisotopic pks. disp.
tmipval = uicontrol('parent',tgA,'Style','text','String',['+/- ',num2str(tmipwindow),...
    ' min(s)'],'Units','normalized','Position',[.28 .655 .03 .05],'FontSize',10);
tmipinc = uicontrol('parent',tgA,'Style','pushbutton','String','+','Units','normalized',...
    'Position',[.28 .71 .03 .05],'Value',1,'FontSize',24,'Callback',@tmipinc_Callback);
tmipred = uicontrol('parent',tgA,'Style','pushbutton','String','-','Units','normalized',...
    'Position',[.28 .58 .03 .05],'Value',1,'FontSize',24,'Callback',@tmipred_Callback);
%% Get necessary parameters and data
scoreintdata_var = displayoptions.scoreintdata;
scoreintdata = scoreintdata_var;
originalsimultaneousfragsetting = scoreintdata_var.sliminput.simultaneousfrag;
numfragmodes = size(displaydataind,2);
current_displaydataind = displaydataind(selectedrow,:);
if isnumeric(current_displaydataind)
    tempcurrent_displaydataind = current_displaydataind;
    current_displaydataind = cell(size(current_displaydataind));
    for j = 1:length(current_displaydataind)
        current_displaydataind{j} = tempcurrent_displaydataind(j);
    end
end
current_numisomers = length(current_displaydataind{1});
currentisomerind = 1;
currentisomerind_sort = 1:current_numisomers;

represent_currentdispind = zeros(current_numisomers,1);  % this is to find a non-zero dispind for each isomer
for j = 1:numfragmodes
    temp_currentdispind = current_displaydataind{j};
    represent_currentdispind(temp_currentdispind > 0) = temp_currentdispind(temp_currentdispind > 0);
end
current_quality3 = usercustom.Quality3(represent_currentdispind,:);

set(prevrow,'visible','off');
if selectedrow > 1
    set(prevrow,'visible','on');
end
if selectedrow == size(displaydataind,1)
    set(nextrow,'visible','off');
end
set(previso,'visible','off');
if current_numisomers == 1
    set(nextiso,'visible','off');
else
    set(nextiso,'visible','on');
end
set(cb_specqual,'Value',current_quality3(currentisomerind,1));
set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
set(cb_subjective,'Value',current_quality3(currentisomerind,3));
set(currentrow,'String',[num2str(selectedrow),'/',num2str(size(displaydataind,1))]);
set(currentiso,'String',['1/',num2str(current_numisomers)]);
%% Initialization is completed

%% Annotate the first set of data
[thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
    structfigsto,currentind,spectrum,theores] = ...
    createtriggeredannotations(msdata,scoredata,...
    current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
    {},{},{},{});
tempcurrentind = currentind(currentind > 0);
current_scoredata= scoredata(tempcurrentind);
plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
    '  CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
    '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
    '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
originalsgp = usercustom.SGPwas{tempcurrentind(1)};
set(editbox_SGPwas,'String',originalsgp);
if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%     set(annotategroup,'String','N/A');
    set(SGPwashelptext,'String','Candidate SGP changed from: ');
    for j = 1:length(annobuttonsto)
        if ~isempty(annobuttonsto{j})
            set(annobuttonsto{j},'String','N/A - SGP changed');
        end
    end
end
    function prevrow_callback(src,event)
        set(nextrow,'visible','on');
        calcing = msgbox('Calculating. Please wait...');
        selectedrow = selectedrow - 1;
        if selectedrow < 1
            selectedrow = 1;
            set(prevrow,'visible','off');
        else
            if selectedrow == 1
                set(prevrow,'visible','off');
            else
                set(nextrow,'visible','on');
            end
            if iscell(displaydataind)
                if length(displaydataind{selectedrow,1}) == 1
                    set(nextiso,'visible','off');
                else
                    set(nextiso,'visible','on');
                end
            elseif isnumeric(displaydataind)
                if sum(sum(displaydataind(selectedrow,:) > 0)) == 1
                    set(nextiso,'visible','off');
                else
                    set(nextiso,'visible','on');
                end
            end
            set(previso,'visible','off');
            currentisomerind = 1;
            current_displaydataind = displaydataind(selectedrow,:);
            if isnumeric(current_displaydataind)
                tempcurrent_displaydataind = current_displaydataind;
                current_displaydataind = cell(size(current_displaydataind));
                for ii = 1:length(current_displaydataind)
                    current_displaydataind{ii} = tempcurrent_displaydataind(ii);
                end
            end
            current_numisomers = length(current_displaydataind{1});
            currentisomerind_sort = 1:current_numisomers;
            set(currentrow,'String',[num2str(selectedrow),'/',num2str(size(displaydataind,1))]);
            set(currentiso,'String',['1/',num2str(current_numisomers)]);
            represent_currentdispind = zeros(current_numisomers,1);  % this is to find a non-zero dispind for each isomer
            for ii = 1:numfragmodes
                temp_currentdispind = current_displaydataind{ii};
                represent_currentdispind(temp_currentdispind > 0) = temp_currentdispind(temp_currentdispind > 0);
            end
            current_quality3 = usercustom.Quality3(represent_currentdispind,:);
            set(cb_specqual,'Value',current_quality3(currentisomerind,1));
            set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
            set(cb_subjective,'Value',current_quality3(currentisomerind,3));
            [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                structfigsto,currentind,spectrum,theores] = ...
                createtriggeredannotations(msdata,scoredata,...
                current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                {},thisfigsto,annobuttonsto,structfigsto);
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            tmipwindow = 0.25;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            originalsgp = usercustom.SGPwas{tempcurrentind(1)};
            set(editbox_SGPwas,'String',originalsgp);
%             set(annotategroup,'String','AnnoAll');
            set(SGPwashelptext,'String','Candidate SGP unchanged: ');
            if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                 set(annotategroup,'String','N/A');
                set(SGPwashelptext,'String','Candidate SGP changed from: ');
                for ii = 1:length(annobuttonsto)
                    if ~isempty(annobuttonsto{ii})
                        set(annobuttonsto{ii},'String','N/A - SGP changed');
                    end
                end
            end
        end
        close(calcing);
    end
    function nextrow_callback(src,event)
        set(prevrow,'visible','on');
        calcing = msgbox('Calculating. Please wait...');
        selectedrow = selectedrow + 1;
        if selectedrow > size(displaydataind,1)
            selectedrow = size(displaydataind,1);
            set(nextrow,'visible','off');
            set(prevrow,'visible','on');
        else
            if selectedrow == size(displaydataind,1)
                set(nextrow,'visible','off');
            else
                set(nextrow,'visible','on');
            end
            if iscell(displaydataind)
                if length(displaydataind{selectedrow,1}) == 1
                    set(nextiso,'visible','off');
                else
                    set(nextiso,'visible','on');
                end
            elseif isnumeric(displaydataind)
                if sum(sum(displaydataind(selectedrow,:) > 0)) == 1
                    set(nextiso,'visible','off');
                else
                    set(nextiso,'visible','on');
                end
            end
            set(previso,'visible','off');
            currentisomerind = 1;
            current_displaydataind = displaydataind(selectedrow,:);
            if isnumeric(current_displaydataind)
                tempcurrent_displaydataind = current_displaydataind;
                current_displaydataind = cell(size(current_displaydataind));
                for ii = 1:length(current_displaydataind)
                    current_displaydataind{ii} = tempcurrent_displaydataind(ii);
                end
            end
            current_numisomers = length(current_displaydataind{1});
            currentisomerind_sort = 1:current_numisomers;
            set(currentrow,'String',[num2str(selectedrow),'/',num2str(size(displaydataind,1))]);
            set(currentiso,'String',['1/',num2str(current_numisomers)]);
            represent_currentdispind = zeros(current_numisomers,1);  % this is to find a non-zero dispind for each isomer
            for ii = 1:numfragmodes
                temp_currentdispind = current_displaydataind{ii};
                represent_currentdispind(temp_currentdispind > 0) = temp_currentdispind(temp_currentdispind > 0);
            end
            current_quality3 = usercustom.Quality3(represent_currentdispind,:);
            set(cb_specqual,'Value',current_quality3(currentisomerind,1));
            set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
            set(cb_subjective,'Value',current_quality3(currentisomerind,3));
            [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                structfigsto,currentind,spectrum,theores] = ...
                createtriggeredannotations(msdata,scoredata,...
                current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                {},thisfigsto,annobuttonsto,structfigsto);
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            tmipwindow = 0.25;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            originalsgp = usercustom.SGPwas{tempcurrentind(1)};
            set(editbox_SGPwas,'String',originalsgp);
%             set(annotategroup,'String','AnnoAll');
            set(SGPwashelptext,'String','Candidate SGP unchanged: ');
            if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                 set(annotategroup,'String','N/A');
                set(SGPwashelptext,'String','Candidate SGP changed from: ');
                for ii = 1:length(annobuttonsto)
                    if ~isempty(annobuttonsto{ii})
                        set(annobuttonsto{ii},'String','N/A - SGP changed');
                    end
                end
            end
        end
        close(calcing);
    end
    function previso_callback(src,event)
        set(nextiso,'visible','on');
        currentisomerind = currentisomerind - 1;
        if currentisomerind < 1
            currentisomerind = 1;
            set(previso,'visible','off');
        else
            if currentisomerind > 1
                set(previso,'visible','on');
            else
                set(previso,'visible','off');
            end
            set(currentiso,'String',[num2str(currentisomerind),'/',num2str(current_numisomers)]);
            set(cb_specqual,'Value',current_quality3(currentisomerind,1));
            set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
            set(cb_subjective,'Value',current_quality3(currentisomerind,3));
            [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                structfigsto,currentind,spectrum,theores] = ...
                createtriggeredannotations(msdata,scoredata,...
                current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                theofragsto,thisfigsto,annobuttonsto,structfigsto);
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            originalsgp = usercustom.SGPwas{tempcurrentind(1)};
            set(editbox_SGPwas,'String',originalsgp);
%             set(annotategroup,'String','AnnoAll');
            set(SGPwashelptext,'String','Candidate SGP unchanged: ');
            if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                 set(annotategroup,'String','N/A');
                set(SGPwashelptext,'String','Candidate SGP changed from: ');
                for i = 1:length(annobuttonsto)
                    if ~isempty(annobuttonsto{i})
                        set(annobuttonsto{i},'String','N/A - SGP changed');
                    end
                end
            end
        end
    end
    function nextiso_callback(src,event)
        set(previso,'visible','on');
        currentisomerind = currentisomerind + 1;
        if currentisomerind > current_numisomers
            currentisomerind = current_numisomers;
            set(nextiso,'visible','off');
        else
            if currentisomerind == current_numisomers
                set(nextiso,'visible','off');
            else
                set(nextiso,'visible','on');
            end
            set(currentiso,'String',[num2str(currentisomerind),'/',num2str(current_numisomers)]);
            set(cb_specqual,'Value',current_quality3(currentisomerind,1));
            set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
            set(cb_subjective,'Value',current_quality3(currentisomerind,3));
            [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                structfigsto,currentind,spectrum,theores] = ...
                createtriggeredannotations(msdata,scoredata,...
                current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                theofragsto,thisfigsto,annobuttonsto,structfigsto);
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            originalsgp = usercustom.SGPwas{tempcurrentind(1)};
            set(editbox_SGPwas,'String',originalsgp);
%             set(annotategroup,'String','AnnoAll');
            set(SGPwashelptext,'String','Candidate SGP unchanged: ');
            if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                 set(annotategroup,'String','N/A');
                set(SGPwashelptext,'String','Candidate SGP changed from: ');
                for i = 1:length(annobuttonsto)
                    if ~isempty(annobuttonsto{i})
                        set(annobuttonsto{i},'String','N/A - SGP changed');
                    end
                end
            end
        end
    end
    function jump_callback(src,event)
        jumptomode = jumptowhat.Value;
        if jumptomode == 3  % isomer
            currentisomerind = floor(str2double(jumptowhere.String));
            if currentisomerind <= 0
                currentisomerind = 1;                
            elseif currentisomerind > current_numisomers
                currentisomerind = current_numisomers;
            end
            set(previso,'visible','off');
            set(nextiso,'visible','off');
            if currentisomerind > 1
                set(previso,'visible','on');
            end
            if currentisomerind < current_numisomers
                set(nextiso,'visible','on');
            end
            set(currentiso,'String',[num2str(currentisomerind),'/',num2str(current_numisomers)]);
            set(cb_specqual,'Value',current_quality3(currentisomerind,1));
            set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
            set(cb_subjective,'Value',current_quality3(currentisomerind,3));
            [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                structfigsto,currentind,spectrum,theores] = ...
                createtriggeredannotations(msdata,scoredata,...
                current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                theofragsto,thisfigsto,annobuttonsto,structfigsto);
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            originalsgp = usercustom.SGPwas{tempcurrentind(1)};
            set(editbox_SGPwas,'String',originalsgp);
%             set(annotategroup,'String','AnnoAll');
            set(SGPwashelptext,'String','Candidate SGP unchanged: ');
            if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                 set(annotategroup,'String','N/A');
                set(SGPwashelptext,'String','Candidate SGP changed from: ');
                for i = 1:length(annobuttonsto)
                    if ~isempty(annobuttonsto{i})
                        set(annobuttonsto{i},'String','N/A - SGP changed');
                    end
                end
            end
        else
            calcing = msgbox('Calculating. Please wait...');
            tgthit = false;
            switch jumptomode
                case 1  % row
                    tgtrowiso = floor(str2num(jumptowhere.String));
                    if length(tgtrowiso) == 1
                        tgtrow = tgtrowiso;
                        tgtiso = 0;
                    elseif length(tgtrowiso) == 2
                        tgtrow = tgtrowiso(1);
                        tgtiso = tgtrowiso(2);
                    end
                    if tgtrow <= size(displaydataind,1)
                        tgthit = true;
                        selectedrow = tgtrow;
                    end
                case 2  % scan
                    tgtscan = floor(str2double(jumptowhere.String));
                    tgtindex = [scoredata.Scan] == tgtscan;
                    tgtiso = 0;
                    if any(tgtindex)
                        tgtfragmode = scoredata(tgtindex).Fragmode;
                        tgtfragmodeind = ismember(scoreintdata_var.colnames,tgtfragmode);
                        tempdpdataind = displaydataind(:,tgtfragmodeind);
                        for i = 1:size(tempdpdataind,1)
                            if ismember(tgtscan,[scoredata(tempdpdataind{i}).Scan])
                                tgthit = true;
                                selectedrow = i;
                                break
                            end
                        end
                    end
            end
            if ~tgthit
                warndlg('No hit found');
            else
                set(prevrow,'visible','off');
                set(nextrow,'visible','off');
                set(previso,'visible','off');
                set(nextiso,'visible','off');
                set(jumptowhere,'String','');
                if selectedrow > 1
                    set(prevrow,'visible','on');
                end
                set(nextrow,'visible','on');
                if selectedrow >= size(displaydataind,1)
                    selectedrow = size(displaydataind,1);
                    set(prevrow,'visible','on');
                    set(nextrow,'visible','off');
                end
                current_displaydataind = displaydataind(selectedrow,:);
                if isnumeric(current_displaydataind)
                    tempcurrent_displaydataind = current_displaydataind;
                    current_displaydataind = cell(size(current_displaydataind));
                    for i = 1:length(current_displaydataind)
                        current_displaydataind{i} = tempcurrent_displaydataind(i);
                    end
                end
                current_numisomers = length(current_displaydataind{1});
                currentisomerind_sort = 1:current_numisomers;
                currentisomerind = 1;
                if tgtiso <= 0 || tgtiso > current_numisomers
                    currentisomerind = 1;
                else
                    currentisomerind = tgtiso;
                end
                set(previso,'visible','off');
                set(nextiso,'visible','off');
                if iscell(displaydataind)
                    if currentisomerind > 1
                        set(previso,'visible','on');
                    end
                    if currentisomerind < current_numisomers
                        set(nextiso,'visible','on');
                    end
                elseif isnumeric(displaydataind)
                    if sum(sum(displaydataind(selectedrow,:) > 0)) == 1
                        set(nextiso,'visible','off');
                    else
                        set(nextiso,'visible','on');
                    end
                end
                set(currentrow,'String',[num2str(selectedrow),'/',num2str(size(displaydataind,1))]);
                set(currentiso,'String',[num2str(currentisomerind),'/',num2str(current_numisomers)]);
                represent_currentdispind = zeros(current_numisomers,1);  % this is to find a non-zero dispind for each isomer
                for i = 1:numfragmodes
                    temp_currentdispind = current_displaydataind{i};
                    represent_currentdispind(temp_currentdispind > 0) = temp_currentdispind(temp_currentdispind > 0);
                end
                current_quality3 = usercustom.Quality3(represent_currentdispind,:);
                set(cb_specqual,'Value',current_quality3(currentisomerind,1));
                set(cb_ms1qual,'Value',current_quality3(currentisomerind,2));
                set(cb_subjective,'Value',current_quality3(currentisomerind,3));
                [thisfigsto,theofragsto,annobuttonsto,currentCombiES,currentDecoyES,...
                    structfigsto,currentind,spectrum,theores] = ...
                    createtriggeredannotations(msdata,scoredata,...
                    current_displaydataind,currentisomerind,scoreintdata_var,tgA,uiscoretable,...
                    {},thisfigsto,annobuttonsto,structfigsto);
                tempcurrentind = currentind(currentind > 0);
                current_scoredata= scoredata(tempcurrentind);
                tmipwindow = 0.25;
                plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
                set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                    '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                    '  SGP = ',current_scoredata(1).SGP,'  Charge = ',num2str(current_scoredata(1).Charge),...
                    '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
                originalsgp = usercustom.SGPwas{tempcurrentind(1)};
                set(editbox_SGPwas,'String',originalsgp);
%                 set(annotategroup,'String','AnnoAll');
                set(SGPwashelptext,'String','Candidate SGP unchanged: ');
                if ~strcmpi(originalsgp,current_scoredata(1).SGP)
%                     set(annotategroup,'String','N/A');
                    set(SGPwashelptext,'String','Candidate SGP changed from: ');
                    for i = 1:length(annobuttonsto)
                        if ~isempty(annobuttonsto{i})
                            set(annobuttonsto{i},'String','N/A - SGP changed');
                        end
                    end
                end
            end
            close(calcing);
        end
    end
%     function annotategroup_callback(src,event)
%         calcing = msgbox('Calculating. Please wait...');
%         allms2tol = displayoptions.scoreintdata.scoreoptions.ms2tol;
%         allms2tolunit = displayoptions.scoreintdata.scoreoptions.ms2tolunit;
%         allfragmodes = upper(displayoptions.scoreintdata.scoreoptions.analyzefragmode);
%         if displayoptions.scoreintdata.scoreoptions.isHCDtrigger
%             dispfragmodes = upper(displayoptions.scoreintdata.colnames);
%         else
%             dispfragmodes = allfragmodes;
%         end
%         for i = 1:length(current_scoredata)
%             genannofig(current_scoredata(i).Scan,spectrum{i},current_scoredata(i),...
%                 allms2tol(ismember(allfragmodes,current_scoredata(i).Fragmode)),...
%                 allms2tolunit{ismember(allfragmodes,current_scoredata(i).Fragmode)},...
%                 msdata,theores{i},theofragsto{currentisomerind,ismember(dispfragmodes,...
%                 current_scoredata(i).Fragmode)},'');
%         end
%         close(calcing);
%     end
    function saveres_callback(src,event)
        [p,f,~] = fileparts(displayoptions.scorefname);
        [outputfilename,outputfolder] = uiputfile('*.mat','Save results to',fullfile(p,[f,'_reviewed']));
        if ~ischar(outputfilename)
            return
        else
            result = scoredata;
            saving_msgbox = msgbox('Saving. Please wait...');
            save(fullfile(outputfolder,outputfilename),'result','scoreintdata','usercustom');
            close(saving_msgbox);
        end
    end
    function cb_specqual_callback(src,event)
        userinput = ~current_quality3(currentisomerind,1);
        current_quality3(:,1) = userinput;
        for ii = 1:numfragmodes
            tempdispdataind = current_displaydataind{ii};
            if tempdispdataind > 0
                usercustom.Quality3(tempdispdataind,1) = userinput;
            end
        end
        set(cb_specqual,'Value',userinput);
    end
    function cb_ms1qual_callback(src,event)
        userinput = ~current_quality3(currentisomerind,2);
        current_quality3(:,2) = userinput;
        for ii = 1:numfragmodes
            tempdispdataind = current_displaydataind{ii};
            if tempdispdataind > 0
                usercustom.Quality3(tempdispdataind,2) = userinput;
            end
        end
        set(cb_ms1qual,'Value',userinput);
    end
    function cb_subjective_callback(src,event)
        userinput = ~current_quality3(currentisomerind,3);
        current_quality3(currentisomerind,3) = userinput;
        for ii = 1:numfragmodes
            tempdispdataind = current_displaydataind{ii}(currentisomerind);
            if tempdispdataind > 0
                usercustom.Quality3(tempdispdataind,3) = userinput;
            end
        end
        set(cb_subjective,'Value',userinput);
        tempcurrent_quality3 = current_quality3(:,3);
        tempcurrentisomerind_sort = currentisomerind_sort;
        tempcurrentisomerind_sort(tempcurrent_quality3 == 0) = ...
            tempcurrentisomerind_sort(tempcurrent_quality3 == 0) + current_numisomers;
        [~,ind] = sort(tempcurrentisomerind_sort);
        currentisomerind_sort = currentisomerind_sort(ind);
        current_quality3 = current_quality3(ind,:);
        for ii = 1:numfragmodes
            current_displaydataind{ii} = current_displaydataind{ii}(ind);
            displaydataind{selectedrow,ii} = displaydataind{selectedrow,ii}(ind);
            theofragsto(:,ii) = theofragsto(ind,ii);
        end
    end
    function rejectall_callback(src,event)
        current_quality3 = false(size(current_quality3));
        for ii = 1:numfragmodes
            tempdispdataind = current_displaydataind{ii};
            if tempdispdataind > 0
                usercustom.Quality3(tempdispdataind,:) = false;
            end
        end
        set(cb_specqual,'Value',false);
        set(cb_ms1qual,'Value',false);
        set(cb_subjective,'Value',false);
    end
    function changecandidate_callback(src,event)
        answer = inputdlg('Input a new SGP sequence to replace current one:','Update Candidate',...
            [1 60],{current_scoredata(1).SGP});
        if ~isempty(answer)
            newsgp = answer{1};
%             set(annotategroup,'String','AnnoAll');
            for ii = 1:length(currentind)
                if currentind(ii) > 0
                    scoredata(currentind(ii)).SGP = newsgp;
                end
            end
            tempcurrentind = currentind(currentind > 0);
            current_scoredata= scoredata(tempcurrentind);
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
            set(tgA,'name',['CombiES = ',num2str(round(currentCombiES,4)),...
                '  ','CombiDecoyES = ',num2str(round(currentDecoyES,4)),...
                '  SGP = ',newsgp,'  Charge = ',num2str(current_scoredata(1).Charge),...
                '  Quant = ',num2str(current_scoredata(1).Quant),'  Protein = ',current_scoredata(1).Protein]);
            for ii = 1:length(currentind)
                if currentind(ii) > 0
                    if strcmpi(usercustom.SGPwas{currentind(ii)},newsgp)
                        set(annobuttonsto{ii},'String','Annotate');
                    else
                        set(annobuttonsto{ii},'String','N/A - SGP changed');
%                         set(annotategroup,'String','N/A');
                    end
                end
            end
        end
    end
    function ssa_callback(src,event)
        spectraAnnotationgui('MSDATA',msdata,'SGP',current_scoredata(1).SGP);
    end
    function suggcandi_callback(src,event)
        suggestisotopicandidates({displayoptions.scoreintdata.sliminput.pepfile},...
            current_scoredata(1).SGP,...
            displayoptions.scoreintdata.scoreoptions.ms1tol,...
            displayoptions.scoreintdata.scoreoptions.ms1tolunit);
    end
    function closereq(src,callbackdata)
        % Close request function
        % to display a question dialog box
        selection = questdlg('Save results before you quit?',...
            'Are you sure you want to quit?',...
            'Save results','Quit Browser','Cancel close','Save results');
        switch selection
            case 'Quit Browser'
                delete(tgA)
            case 'Save results'
                [p,f,~] = fileparts(displayoptions.scorefname);
                [outputfilename,outputfolder] = uiputfile('*.mat','Save results to',fullfile(p,[f,'_reviewed']));
                if ~ischar(outputfilename)
                    return
                else
                    result = scoredata;
                    saving_msgbox = msgbox('Saving. Please wait...');
                    save(fullfile(outputfolder,outputfilename),'result','scoreintdata','usercustom');
                    close(saving_msgbox);
                    delete(tgA)
                end
            case 'Cancel close'
                return
        end
    end
    function tmipinc_Callback(src,event)
        if tmipwindow == 0.25 || tmipwindow == 0.5
            tmipwindow = tmipwindow * 2;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        elseif tmipwindow == 0
            tmipwindow = 0.25;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        elseif tmipwindow == 10
            % INTENTIONALLY LET BLANK
        else
            tmipwindow = tmipwindow + 1;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        end
    end
    function tmipred_Callback(src,event)
        if tmipwindow == 1 || tmipwindow == 0.5
            tmipwindow = tmipwindow / 2;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        elseif tmipwindow == 0.25
            tmipwindow = 0;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        elseif tmipwindow == 0
            tmipwindow = 0;
        else
            tmipwindow = tmipwindow - 1;
            plotmonoisodist(msdata,current_scoredata,MS1dist,tmipval,tmipwindow);
        end
    end
    function simufrag_callback(src,event)
        if simufragstat
            simufragstat = false;
            set(txt_simultaneousfrag,'String','Default simultaneous frag. setting');
            scoreintdata_var.scoreoptions.simultaneousfrag = ...
                scoreintdata.scoreoptions.simultaneousfrag;
            scoreintdata_var.sliminput.simultaneousfrag = ...
                scoreintdata.sliminput.simultaneousfrag;
            scoreintdata_var.scoreoptions.maxstublen = ...
                scoreintdata.scoreoptions.maxstublen;
            scoreintdata_var.sliminput.maxstublen = ...
                scoreintdata.sliminput.maxstublen;
        else
            simufragstat = true;
            residuesize = inputdlg('Enter glycan residue size limit','Pep/Gly simultaneous fragmentation',...
                [1 35],{'5'});
            residuesize = floor(abs(str2double(residuesize{1})));
            set(txt_simultaneousfrag,'String',['Enable simultaneous frag., residue size <= ',num2str(residuesize)]);
            scoreintdata_var.scoreoptions.simultaneousfrag = ...
                true(size(scoreintdata.scoreoptions.simultaneousfrag));
            scoreintdata_var.sliminput.simultaneousfrag = ...
                true(size(scoreintdata.sliminput.simultaneousfrag));
            scoreintdata_var.scoreoptions.maxstublen = ...
                residuesize;
            scoreintdata_var.sliminput.maxstublen = ...
                residuesize;
        end
    end
    function keyboardaction(src,event)
        switch event.Key
            case {'s','S'}
                cb_specqual_callback;
            case {'m','M'}
                cb_ms1qual_callback;
            case {'u','U'}
                cb_subjective_callback;
            case {'r','R'}
                rejectall_callback;
            case 'uparrow'
                prevrow_callback;
            case 'downarrow'
                nextrow_callback;
            case 'leftarrow'
                if current_numisomers > 1
                    previso_callback;
                end
            case 'rightarrow'
                if current_numisomers > 1
                    nextiso_callback;
                end
            case 'a'
%                 annotategroup_callback;
            case {'1','2','3','4','5','6','7','8','9'}
                allms2tol = displayoptions.scoreintdata.scoreoptions.ms2tol;
                allms2tolunit = displayoptions.scoreintdata.scoreoptions.ms2tolunit;
                allfragmodes = upper(displayoptions.scoreintdata.scoreoptions.analyzefragmode);
                if displayoptions.scoreintdata.scoreoptions.isHCDtrigger
                    dispfragmodes = upper(displayoptions.scoreintdata.colnames);
                else
                    dispfragmodes = allfragmodes;
                end
                whichspectrum = str2double(event.Key);
                genannofig(current_scoredata(whichspectrum).Scan,spectrum{whichspectrum},...
                    current_scoredata(whichspectrum),...
                    allms2tol(strcmpi(allfragmodes,current_scoredata(whichspectrum).Fragmode)),...
                    allms2tolunit{strcmpi(allfragmodes,current_scoredata(whichspectrum).Fragmode)},...
                    msdata,theores{whichspectrum},theofragsto{currentisomerind,strcmpi(dispfragmodes,...
                    current_scoredata(whichspectrum).Fragmode)},'');
            case 'pageup'
                tmipinc_Callback;
            case 'pagedown'
                tmipred_Callback;
            otherwise
                % INTENTIONALLY LEFT BLANK
        end
    end
end