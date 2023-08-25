function varargout = browsegui(varargin)
% BROWSEGUI: Start a GUI for browsing scoring result in MAT format
%
% See also: GLYCOPATGUI, DIGESTGUI, SCOREGUI, FRAGGUI, SPECTRAANNOTATIONGUI

% Author: Gang Liu, Kai Cheng
% Date Lastly Updated: 12/20/19

% Edit the above text to modify the response to help browsegui

% Last Modified by GUIDE v2.5 15-Aug-2023 10:32:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @browsegui_OpeningFcn, ...
    'gui_OutputFcn',  @browsegui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% USER EDITABLE REGION START

%% BUTTON - LOAD RESULT FILE
function pushbutton_loadscoreresult_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadscoreresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_clear_Callback(hObject, eventdata, handles)
[dirallfilename, pathname,~] = uigetfile({'*.mat'},'Pick a file storing scoring result');
if(isequal(dirallfilename,0) || isequal(pathname,0))
    errordlg('File NOT Selected','File NOT Selected');
    return
end
handles.displayoptions.scorefname = fullfile(pathname,dirallfilename);
nowloadingmsgbox = msgbox('Now Loading...');
% read header from file
[~,~,ext] = fileparts(handles.displayoptions.scorefname);
if strcmpi(ext,'.mat')
    resultstruct = load(handles.displayoptions.scorefname);
    scoredata = resultstruct.result;
    scoreintdata = resultstruct.scoreintdata;
    if isfield(handles.displayoptions,'scorefname')
        if contains(handles.displayoptions.scorefname,'_O_glycopep_1_1.mat')
            load('allfragments_O12.mat','allfragments');
        else
            load('allfragments_N1008.mat','allfragments');
        end
    end
    % if ~isfield(scoreintdata.sliminput,'allfragments_path')|| ...
    %         ~exist(scoreintdata.sliminput.allfragments_path,'file')
    %     [file,path] = uigetfile('.mat','Select PTM Fragment File');
    %     allfragmentspath = fullfile(path,file);
    % else
    %     allfragmentspath = scoreintdata.sliminput.allfragments_path;
    % end
    % load(allfragmentspath,'allfragments');
    allfragments_fldnms = fieldnames(allfragments);
    for ii = 1:length(allfragments_fldnms)
        scoreintdata.(allfragments_fldnms{ii}) = allfragments.(allfragments_fldnms{ii});
    end
    allptmseq = allfragments.ptmseq;
    for ii = 1:length(scoredata)
        tempprotid = str2num(scoredata(ii).ProteinID);
        tempsgp = scoredata(ii).SGP;
        [~,g,m] = breakGlyPep(tempsgp);
        tempptmseq = {};
        tempptmpos = [];
        if ~isempty(g)
            tempptmseq = [tempptmseq,g.struct];
            tempptmpos = [tempptmpos,g.pos];
        end
        if ~isempty(m)
            tempptmseq = [tempptmseq,m.struct];
            tempptmpos = [tempptmpos,m.pos];
        end
        [~,ind] = sort(tempptmpos);
        tempptmseq = tempptmseq(ind);
        for jj = 1:length(tempptmseq)
            tempprotid(jj + 2) = find(strcmpi(allptmseq,tempptmseq{jj}));
        end
        scoredata(ii).ProteinID = num2str(tempprotid);
    end




    handles.displayind = [];
    handles.displayoptions.scoreintdata = scoreintdata;
    
    scoredata_examine = cell(size(scoredata));
    for ii = 1:length(scoredata)
        scoredata_examine{ii} = [num2str(scoredata(ii).Scan),' ',scoredata(ii).ProteinID];
    end
    [~,uniind] = unique(scoredata_examine);
    scoredata = scoredata(uniind);    
    
    handles.displayoptions.exptdatafname = scoreintdata.sliminput.exptdata;
%    handles.displayoptions.exptdatafname = [scoreintdata.sliminput.exptdata,'.mat'];  % SN edit
    set(handles.edit_exptdata,'String',handles.displayoptions.exptdatafname);   % SN edit
    handles.usercustom = [];  % a wrapper for user subjective selection & auxiliary data
    if isfield(resultstruct,'usercustom')
        handles.usercustom = resultstruct.usercustom;
    else  % this happens when loading data for the first time, "build usercustom here"
        % In the future when adding new fields, initialize here.
        handles.usercustom.Quality3 = true(length(scoredata),3);
        tempscoresgps = {scoredata.SGP};
        handles.usercustom.SGPwas = tempscoresgps(:);
    end
    handles.scoredata = scoredata;
    set(handles.edit_loadfilestatus,'String',fullfile(pathname,dirallfilename));
    set(handles.edit_dispheader,'String',compileheadertext(scoreintdata));
else
    errordlg('Unsupported File Type.')
    return
end
if isfield(handles,'msdata')
    handles = rmfield(handles,'msdata');
end
close(nowloadingmsgbox);
%combineisomer = questdlg('Combine isomers?','','Yes','No','Yes');
combineisomer = 'Yes'
switch combineisomer
    case 'Yes'
        handles.displayoptions.combineisomer = true;
    case 'No'
        handles.displayoptions.combineisomer = false;
end
calcing = msgbox('Please wait...');
[displaydata,displaydataind,displaycolnames,displaycolformat] = builddisptable(...
    scoredata,[],handles.displayoptions,handles.usercustom);
handles.displaydataind = displaydataind;
handles.displaydataind_original = displaydataind;
set(handles.uitable_scoredisplay,'data',displaydata,'ColumnName',...
    displaycolnames,'ColumnFormat',displaycolformat);
% Previous section has the table ready for display
% Need to adjust filter settings according to experiment data
unifragmode = unique({scoredata.Fragmode});
if scoreintdata.scoreoptions.isHCDtrigger
    fdrfragmodes = cell(length(unifragmode)+3,1);
    applytofragmodes = cell(length(unifragmode)+3,1);
    fdrfragmodes{1} = 'CombiES';
    fdrfragmodes{2} = 'Every hit';
    fdrfragmodes{3} = 'Each frag. individually';
    applytofragmodes{1} = 'CombiES';
    applytofragmodes{2} = 'Everyone';
    applytofragmodes{3} = 'Each frag. ind.';
    fdrfragmodes(4:end) = unifragmode;
    applytofragmodes(4:end) = unifragmode;
else
    fdrfragmodes = cell(length(unifragmode) + 2,1);
    applytofragmodes = cell(length(unifragmode) + 2,1);
    fdrfragmodes{1} = 'Every hit';
    fdrfragmodes{2} = 'Each frag. individually';
    applytofragmodes{1} = 'Everyone';
    applytofragmodes{2} = 'Each frag. ind.';
    fdrfragmodes(3:end) = unifragmode;
    applytofragmodes(3:end) = unifragmode;
end
set(handles.fdr_fragmode,'String',fdrfragmodes);
set(handles.select_fragmode,'String',applytofragmodes);
close(calcing);
guidata(hObject,handles);


%% BUTTON - RELOAD .MAT FILE
function pb_loaddatamat_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loaddatamat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'msdata')
    handles = rmfield(handles,'msdata');
end
handles.displayoptions.exptdatafname = '';
[dirallfilename, pathname,~] = uigetfile({'*.mat'},'Pick the file storing experiment data');
handles.displayoptions.exptdatafname = fullfile(pathname,dirallfilename);
guidata(hObject,handles);


%% CELL SELECTION
function uitable_scoredisplay_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_scoredisplay (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.Indices)
    if isfield(handles,'msdata')
        msdata = handles.msdata;
    else
        msdata = [];
    end
    selectedrow = eventdata.Indices(1);
    if isempty(msdata)
        datafname = handles.displayoptions.exptdatafname;
        if ~endsWith(datafname,'.mat')
            datafname=[datafname,'.mat']; %SN edit
        end
        if ~exist(datafname,'file')
                warndlg('Data file not present','File not exist')
            return
        end
        fileloading = msgbox("Now loading...","File Loading");
        msdatastruct = load(datafname);
        fldnm = fieldnames(msdatastruct);
        msdata = msdatastruct.(fldnm{1});
        fileloading = msgbox('Completed.','File Loading','replace');
        load('default_iso_dist.mat','defaultisodist');
        msdata.defaultisodist = defaultisodist
        if ~isfield(msdata,'defaultisodist')
            [file,path] = uigetfile('.mat','Select Isotope Model File');
            defaultisodistpath = fullfile(path,file);
            load(defaultisodistpath,'defaultisodist');
            msdata.defaultisodist = defaultisodist;
        end
        pause(0.1);
        clear fileloading;
        handles.msdata = msdata;
    end
    scoredata = handles.scoredata;
    displaydataind = handles.displaydataind;
    displayoptions = handles.displayoptions;
    usercustom = handles.usercustom;
    annotatespectra(msdata,scoredata,displaydataind,displayoptions,usercustom,selectedrow);
end
guidata(hObject,handles);


%% BUTTON - CLEAR ALL
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_loadfilestatus,'String','');
set(handles.edit_exptdata,'String','');
set(handles.edit_dispheader,'String','');
set(handles.uitable_scoredisplay,'data',[],'ColumnName','');
set(handles.fdr_rate,'Value',1);
set(handles.fdr_fragmode,'String',{'All Modes'},'Value',1);
set(handles.escutoff,'String','');
set(handles.checkbox_selectspectrabyescore,'Value',0);
set(handles.checkbox_precionselection,'Value',0);
set(handles.checkbox_spectraquality,'Value',0);
set(handles.checkbox_userdecision,'Value',0);
set(handles.enscorethreshold,'String','');
set(handles.select_fragmode,'String',{'All Modes'},'Value',1);
if isfield(handles,'scoredata')
    handles = rmfield(handles,'scoredata');
end
if isfield(handles,'displayind')
    handles = rmfield(handles,'displayind');
end
if isfield(handles,'displayoptions')
    handles = rmfield(handles,'displayoptions');
end
if isfield(handles,'usercustom')
    handles = rmfield(handles,'usercustom');
end
clear global
guidata(hObject,handles);



%% BUTTON - EXPORT CURRENT TABLE TO CSV FILE
function export_csv_Callback(hObject, eventdata, handles)
% hObject    handle to export_csv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,p] = uiputfile('*.xlsx','Save Currently Displayed Table to:');
exportdatafilename = fullfile(p,f);
 scoreCSVwrite(exportdatafilename,handles);


%% BUTTON - RESET TO ORIGINAL
function resetoriginal_Callback(hObject, eventdata, handles)
% hObject    handle to resetoriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable_scoredisplay,'data',[]);
calcing = msgbox('Please wait...');
[displaydata,displaydataind,displaycolnames,displaycolformat] = builddisptable(...
    handles.scoredata,handles.displaydataind_original,handles.displayoptions,handles.usercustom);
handles.displaydataind = displaydataind;
set(handles.uitable_scoredisplay,'data',displaydata,'ColumnName',...
    displaycolnames,'ColumnFormat',displaycolformat);
close(calcing);
guidata(hObject,handles);


%% BUTTON - CALCULATE ES CUTOFF
function fdrfilter_Callback(hObject, eventdata, handles)
% hObject    handle to fdrfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set parameters to default
set(handles.escutoff,'String','');
set(handles.enscorethreshold,'String','');
set(handles.checkbox_selectspectrabyescore,'Value',0);
set(handles.select_fragmode,'Value',1);

% Retrieve user FDR level selection
allfdrlvl = cellfun(@(x) str2double(x(1:end-1)), get(handles.fdr_rate,'String'))/100;
fdrlvlind = get(handles.fdr_rate,'Value');
fdrfiltervalue = allfdrlvl(fdrlvlind);

allfdrfragmodes = get(handles.fdr_fragmode,'String');
fdrfragmodeind = get(handles.fdr_fragmode,'Value');
fdrfragmode = allfdrfragmodes{fdrfragmodeind};

escutoff = fdrfilter(handles.scoredata,handles.displayoptions,...
    handles.displaydataind,fdrfiltervalue,fdrfragmode);
if length(escutoff) > 1
    fragmethods = handles.displayoptions.scoreintdata.scoreoptions.analyzefragmode;
    escutoffstring = '';
    for i = 1:length(escutoff)
        escutoffstring = [escutoffstring,fragmethods{i},':',num2str(escutoff(i)),','];
    end
    escutoffstring = escutoffstring(1:end-1);
else
    escutoffstring = num2str(escutoff);
end
set(handles.escutoff,'String',escutoffstring);
set(handles.enscorethreshold,'String',escutoffstring);
set(handles.checkbox_selectspectrabyescore,'Value',1);
set(handles.checkbox_spectraquality,'Value',1);
set(handles.checkbox_precionselection,'Value',1);
set(handles.checkbox_userdecision,'Value',1);
set(handles.select_fragmode,'Value',fdrfragmodeind);
guidata(hObject,handles);


%% BUTTON - EXPORT ANNOTATED FIGURE
function export_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to export_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createannopdf(handles.msdata,handles.scoredata,handles.displaydataind,...
    handles.displayoptions);



%% BUTTON - APPLY ES CUTOFF
% --- Executes on button press in pb_applyescutoff.
function pb_applyescutoff_Callback(hObject, eventdata, handles)
% hObject    handle to pb_applyescutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable_scoredisplay,'data',[]);
pause(0.001)
calcing = msgbox('Please wait...');
allfilterfragmodes = get(handles.select_fragmode,'String');
filterfragmodeind = get(handles.select_fragmode,'Value');
filterbyfragmode = allfilterfragmodes{filterfragmodeind};
filterbyglypeponly = false;
filterbyes = logical(get(handles.checkbox_selectspectrabyescore,'Value'));
filterbyprecion = logical(get(handles.checkbox_precionselection,'Value'));
filterbyuserdec = logical(get(handles.checkbox_userdecision,'Value'));
filterbyspecqual = logical(get(handles.checkbox_spectraquality,'Value'));
filterbyesstring = get(handles.enscorethreshold,'String');
if ~any(strfind(filterbyesstring,','))
    filterbyescutoff = str2num(filterbyesstring);
    filterbyescutoff = reshape(filterbyescutoff,1,[]);
    if filterbyes && ...
            (any(isnan(filterbyescutoff)) || any(isempty(filterbyescutoff)) || any(filterbyescutoff > 1))
        errordlg('Unsupported enseble score cut-off value',...
            'Unsupported enseble score cut-off value');
        return
    end
else
    filterbyescutoff = cellfun(@str2double,strsplit(filterbyesstring,{',',':'}));
    if mod(filterbyescutoff,2)
        % This means there are multiple es cutoff inputs and Name-Value are
        % not matched properly
        errordlg('Incorrect enseble score cut-off value',...
            'Incorrect enseble score cut-off value');
        return        
    else
        filterbyescutoff = filterbyescutoff(2:2:end);
    end
end
[displaydata,displaydataind] = rebuilddisptable(handles.scoredata,...
    handles.displaydataind_original,handles.displayoptions,handles.usercustom,...
    {filterbyfragmode,filterbyglypeponly,filterbyes,filterbyescutoff,...
    filterbyspecqual,filterbyprecion,filterbyuserdec});
if isempty(displaydata)
    warndlg('No hits left!')
end
set(handles.uitable_scoredisplay,'data',displaydata);
handles.displaydataind = displaydataind;
close(calcing);
guidata(hObject,handles);


%% FUNCTION - COMPILE SUMMARY TEXT
function headertext = compileheadertext(scoreintdata)
fragmodes = scoreintdata.scoreoptions.analyzefragmode;
ms2tols = scoreintdata.scoreoptions.ms2tol;
ms2tolunits = scoreintdata.scoreoptions.ms2tolunit;
fragmodetext = '';
for i = 1:length(ms2tols)
    fragmodetext = [fragmodetext,', ',fragmodes{i},': ',num2str(ms2tols(i)),' ',ms2tolunits{i}];    
end
fragmodetext = fragmodetext(3:end);
headertext = {['Peptide file name: ',scoreintdata.sliminput.pepfile];...
    ['Expt. data file name: ',scoreintdata.sliminput.exptdata];...
    ['MS1 tolerence & unit: ',num2str(scoreintdata.scoreoptions.ms1tol),' ',...
    scoreintdata.scoreoptions.ms1tolunit];...
    'Fragmentation mode (Mode, MS2 tolerence & unit):';fragmodetext;...
    ['Cross Correlation maxLag: ',num2str(scoreintdata.sliminput.maxlag)];...
    ['Spectra Cleanup param (CutOffMed and FracMax)= ',...
    num2str(reshape(scoreintdata.sliminput.cutoffmed,1,[])),', ',num2str(reshape(scoreintdata.sliminput.fracmax,1,[]))];...
    ['Select Peak = ',num2str(scoreintdata.sliminput.selectpeak')]};

%% USER EDITABLE REGION END


%% BELOW ARE MATLAB SYSTEM FUNCTION - DO NOT EDIT!

function browsegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to browsegui (see VARARGIN)

% Choose default command line output for browsegui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function varargout = browsegui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function enscorethreshold_Callback(hObject, eventdata, handles)
% hObject    handle to enscorethreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enscorethreshold as text
%        str2double(get(hObject,'String')) returns contents of enscorethreshold as a double


function edit_loadfilestatus_Callback(hObject, eventdata, handles)
% hObject    handle to edit_loadfilestatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_loadfilestatus as text
%        str2double(get(hObject,'String')) returns contents of edit_loadfilestatus as a double


% --- Executes during object creation, after setting all properties.
function edit_loadfilestatus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_loadfilestatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dispheader_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dispheader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dispheader as text
%        str2double(get(hObject,'String')) returns contents of edit_dispheader as a double


% --- Executes during object creation, after setting all properties.
function edit_dispheader_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dispheader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable_scoredisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_scoredisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when annotation is resized.
function annotation_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_loadscoreresult.
function pushbutton_loadscoreresult_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadscoreresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in uitable_scoredisplay.
function uitable_scoredisplay_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_scoredisplay (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function enscorethreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enscorethreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in select_fragmode.
function select_fragmode_Callback(hObject, eventdata, handles)
% hObject    handle to select_fragmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if ~handles.scoreintdata.sliminput.dohcdtrigger
%
% end
% Hints: contents = cellstr(get(hObject,'String')) returns select_fragmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_fragmode


% --- Executes during object creation, after setting all properties.
function select_fragmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_fragmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fdr_rate.
function fdr_rate_Callback(hObject, eventdata, handles)
% hObject    handle to fdr_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function fdr_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fdr_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fdr_fragmode.
function fdr_fragmode_Callback(hObject, eventdata, handles)
% hObject    handle to fdr_fragmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fdr_fragmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fdr_fragmode


% --- Executes during object creation, after setting all properties.
function fdr_fragmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fdr_fragmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fdr_peptidetype.
function fdr_peptidetype_Callback(hObject, eventdata, handles)
% hObject    handle to fdr_peptidetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fdr_peptidetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fdr_peptidetype


% --- Executes during object creation, after setting all properties.
function fdr_peptidetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fdr_peptidetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function escutoff_Callback(hObject, eventdata, handles)
% hObject    handle to escutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of escutoff as text
%        str2double(get(hObject,'String')) returns contents of escutoff as a double


% --- Executes during object creation, after setting all properties.
function escutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to escutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fdr_peptype.
function fdr_peptype_Callback(hObject, eventdata, handles)
% hObject    handle to fdr_peptype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function fdr_peptype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fdr_peptype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on uitable_scoredisplay and none of its controls.
function uitable_scoredisplay_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable_scoredisplay (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton_loadscoreresult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadscoreresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close main.
function main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% check if there's new quant data remains, ask if user wants to save
% Hint: delete(hObject) closes the figure
delete(hObject);



% --- Executes on button press in checkbox_spectraquality.
function checkbox_spectraquality_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_spectraquality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_spectraquality


% --- Executes on button press in checkbox_precionselection.
function checkbox_precionselection_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_precionselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_precionselection


% --- Executes on button press in checkbox_userdecision.
function checkbox_userdecision_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_userdecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_userdecision

%% CHECKBOX - SHOW GLYPEP ONLY
function checkbox_glycopeponly_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_glycopeponly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% CHECKBOX - SET ENSEMBLE SCORE CUTOFF
function checkbox_selectspectrabyescore_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_selectspectrabyescore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in edit_exptdata.
function edit_exptdata_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exptdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns edit_exptdata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_exptdata


% --- Executes during object creation, after setting all properties.
function edit_exptdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exptdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
