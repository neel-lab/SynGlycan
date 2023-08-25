function statusreporthandles = scoreAllSpectra_test(msdata,input)
% SCOREALLSPECTRA: Main function for spectra matching
%
% Syntax:
% statusreporthandles = scoreAllSpectra(msdata,input)
%
% Input:
% msdata: structure. Mass spectrometry experiment data.
% In order for the program to function correctly, msdata must contain the
%     following fields:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     scannum                            n x 1 double                        Scan number of spectra
%     mslvl                                    n x 1 double                       MS level of scan
%     retime                                 n x 1 double                       Retention time of scan
%     fragmode                           n x 1 cell array of              Fragmentation mode of scan
%                                                       strings
%     precursormz                      n x 1 double                      Monoisotopic m/z of precursor ion
%     allprecursormz                 n x 1 cell array of              Precursor m/z information, the 3 values are:
%                                                       1 x 3 double                      [selected ion,isolation window target, monosiotopic m/z].
%     charge                                 n x 1 double                       Charge of precursor ion
%     precursorScanNum         n x 1 double                       Precursor  scan number of MS2 scans, value is -1 for MS1 scans
%     spectra                                n x 1 cell array of              Spectra, 1st column is m/z, 2nd column is intensity
%                                                       m x 2 double
%     parsedAUC                         n x 1 double                       AUC calculated from precursor m/z.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
% input: structure. All other information necessary for scoring.
% The fields are:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     pepfile                                String                                   Full path of digested glycopeptide file.
%     exptdata                            String                                   Full path of processed experiment data file (msdata).
%     fragmode                          n x 1 cell array of              MS2 fragmentation modes that user wants to analyze.
%                                                    strings
%     ms1tol                                Double                                 Tolerence of MS1 matching (value).
%     ms1tolunit                         String                                   Tolerence of MS1 matching (unit).
%     ms2tol                                n x 1 double                       Tolerence of MS2 matching corresponding to
%                                                                                                      each fragmentation mode (value).
%     ms2tolunit                         n x 1 cell array of              Tolerence of MS2 matching corresponding to
%                                                     strings                                    each fragmentation mode (unit).
%     minmaxmz                         1 x 2 double                       Lower and upper limit of m/z values of peaks in spectrum to be
%                                                                                                     matched.
%     outputdir                           String                                    Path where the result files will be saved.
%     outputfname                    String                                   Template of result file name. Actual file names will be extended
%                                                                                                    by the index number of proteins for distinguishing.
%     maxlag                                Double                                 Maximum lag for cross correlation analysis.
%     cutoffmed                         n x 1 Double                       Denoising factor. See SPECTRADENOISING for detail.
%     fracmax                              n x 1 Double                       Denoising factor. See SPECTRADENOISING for detail.
%     fragnum                             n x 3 double                        Maximum number of fragmentations allowed on each part of
%                                                                                                     the glycopeptide. The 3 parts are:
%                                                                                                     [Peptide, Glycan, Non-glycan PTM].
%     fragiontyp                          n x 5 logical                         Peptide fragment ion types to be considered during analysis.
%                                                                                                      The 5 ion types are: [b-, c-, y-, z-, i-].
%     selectpeak                          n x 1 double                        During each spectrum analysis, the existence of these specific
%                                                                                                      m/z values will be checked.
%     monosacislabile                n x 1 logical                         If is true, when creating theoretical fragment ions, if the ion
%                                                                                                       contains Fuc or Neu5Ac, these ions will be further fragmented
%                                                                                                       where these 2 monosaccharides will detach.
%     handles                                Struct                                   The GUI handle of scoregui. Provides access to components
%                                                                                                       of GUI.
%     numprotperworker          Double                                 The maximum number of proteins each worker will handle
%                                                                                                       each time.
%     dohcdtrigger                       Logical                                 Whether this experiment is HCD-triggered. If true, an additional
%                                                                                                       experiment file containing information about association
%                                                                                                       between scans and fragmentation methods is required.
%     maxstublen                         Double                                 The maximum length of glycan residue allowed on peptide
%                                                                                                        fragments.
%     simultaneousfrag              n x 1 logical                        Whether peptide/glycan simultaneous fragmentation will be
%                                                                                                        considered, if both number of cleavages are greater than zero.
%     addoxoniumion                 n x 1 logical.                       If true, when specific fragment (Neu5Ac, Hex, HexNAc) appears,
%                                                                                                        these ions will be further fragmented where small, neutral
%                                                                                                        molecules (water, C2H4O2, CH6O3) will detach.
%     doparacomp                       Double                                The number of workers (CPU cores) to be used.
%     presearch                            Logical                                 Program will use a built-in small glycan library to perform a
%                                                                                                        smaller scale search. Only the proteins identified in this search
%                                                                                                        will be adopted for the full scale search.
%     presearchmode                  Double                               Algorithm to use during presearch. 
%                                                                                                      1 means program will search for a pre-determined set of glycans.
%                                                                                                          Glycoproteins that are confidently identified will be searched
%                                                                                                          for.
%                                                                                                      2 means the program will search for all glycopeptides, for a PSM
%                                                                                                          be included in subsequent steps, Top10 value must >= 3.
%                                                                                                      3 is similar to 2. Criteria changed to at least >= 2 peptide bond 
%                                                                                                          fragmentations be found.
%     presearchtopprotnum     Double                                If presearch is to be performed, how many proteins should be
%                                                                                                       included in the main search.
%     resumeprevsearch            Logical                                 Whether or not check existing files in output folder to decide
%                                                                                                      where the search be started. This is useful for restarting a
%                                                                                                      search that was previously interrupted. By default this is set to
%                                                                                                      TRUE.
%     analyzetgt                            Double                                Type of candidate to be analyzed. 1 is analyze glycopeptide only,
%                                                                                                      2 is analyze bothe glycopeptide and peptide, 3 is analyze
%                                                                                                      peptide only.
%     analyglytyp                          String                                  Type of glycans to analyze. Currently supports "N" and "O". Does
%                                                                                                      not matter when analyzing peptides.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%
% Output:
% statusreporthandles: A copy of input.handles, this is the channel between
%     scoregui and this program.
% The scoring result will be saved as .mat files. The full name of these
%     files is as follows: outputdir\outputfname_ID of the first protein
%     searched in this batch_ID of the last protein searched in this batch.mat
% For example, in file "D:\result\Analysis_1_100.mat", "D:\result\" is
%     where all result files will be saved, defined in "input.outputdir".
%     "Analysis" is the name template for all result files, defined in
%     "input.outputfname". "1" and "100" means the analysis result of
%     protein No. 1 through No. 100 is stored in this file. These 2 numbers
%     are calculated from "input.numprotperworker" and "input.doparacomp".
%
% Note:
% See children function SEARCHSPECTRA for more details about program.
%
% Example:
% N/A. Please run the program with test files provided and set break
%     points.
% To test function: load peptide file
%     '"Test Data Folder"\Scoring\DIGESTED.txt'.
% Then load experiment data file
%     '"Test Data Folder"\Scoring\TESTDATA.mat'.
% Click "Score" to run.
%
% Children function:
% N/A
%
% See Also:
% SEARCHSPECTRA DIGESTFILEANALY ANALYZEEXISTINGFILES ANALYZEPRESEARCHDATA
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

tic  % Start timer

%% Step 1: Unpack inputs and set up parameters
pepfile = input.pepfile;
exptdata = input.exptdata;
scoreoptions.analyzefragmode = input.fragmode;
scoreoptions.ms1tol = input.ms1tol;
scoreoptions.ms1tolunit = input.ms1tolunit;
scoreoptions.ms2tol = input.ms2tol;
scoreoptions.ms2tolunit = input.ms2tolunit;
scoreoptions.minmaxmz = input.minmaxmz;
scoreoptions.maxlag = input.maxlag;
scoreoptions.cutoffmed = input.cutoffmed;
scoreoptions.fracmax = input.fracmax;
scoreoptions.fragnum = input.fragnum;
scoreoptions.selectpeak = input.selectpeak;
scoreoptions.assemblegpmode = 1;
scoreoptions.numprotperworker = input.numprotperworker;
scoreoptions.resumeprevsearch = input.resumeprevsearch;
scoreoptions.analyzeglycan = input.analyzeglycan;
scoreoptions.maxOglyonpep = input.maxOglyonpep;

% Special options
scoreoptions.doparallel = false;
scoreoptions.monosacislabile = input.monosacislabile;
scoreoptions.simultaneousfrag = input.simultaneousfrag;
scoreoptions.addoxoniumion = input.addoxoniumion;
scoreoptions.isHCDtrigger = input.dohcdtrigger;
scoreoptions.maxstublen = input.maxstublen;
scoreoptions.presearch = input.presearch;

% HCD triggered experiment needs additional input
% Trigger data ("SASSO", "colnames") are stored in a separate file, need to be loaded before use
% Filename for trigger data must be expt. data's file name with "_trigger"
% e.g. expt. data file name is "X.mat", trigger data file name has to be
%     "X_trigger.mat", otherwise file will not be read.
triggerdata = [];
if input.dohcdtrigger  % HCD trigger - supplemental activation experiments
    [path,filename,~] = fileparts(exptdata);
    sassofilename = fullfile(path,[filename,'_trigger.mat']);
    triggerdata = load(sassofilename);
end

% GUI handle. It is not essential for the program to operate, can be left
%     blank.
statusreporthandles = [];
statusstr = {};
if isfield(input,'handles')
    scoreguihandles = input.handles;
    statusreporthandles = guidata(scoreguihandles.edit_statusreport);
    statusstr = get(statusreporthandles.edit_statusreport,'string');
end

% Set up peptide fragmentation parameters - inputs are logical values. They
%     are converted to text formats so that downstream functions can get
%     this info directly.
% The fragmentation parameter here is whether or not generating b-, c-, y-,
%     z- and i- ions.
defaultfragiontyp = 'bcyzi';
fragiontyp = input.fragiontyp;
userpepiontyp = cell(size(fragiontyp,1),1);
for ii = 1:length(userpepiontyp)
    userpepiontyp{ii} = defaultfragiontyp(logical(fragiontyp(ii,:)));
end
scoreoptions.userpepiontyp = userpepiontyp;

%% Step 2: Read digested glycoepeptide candidates
% These information include:
%     Protein names;
%     Digested peptide sequences;
%     All PTMs;
%     Fragmentation parameters;
%     Other miscellaneous information: missed cleavage, min. max. length of
%         peptide, number of PTMs allowed, etc.
if ~isempty(statusreporthandles)
    statusstr = [statusstr;'Analyzing digestion file...'];
    set(statusreporthandles.edit_statusreport, 'String', statusstr);
    pause(0.0001);  % leave reaction time for GUI data update
end
allpepdata = digestfileanaly(pepfile,1);  % version 0: GPAT 1; version 1: GPAT 2; version 3: GPAT2+
%% Apply target candidate type setting - Glycopeptide only / Glycopeptide + peptide / Peptide only
% Digested proteins that do not contain target candidates will be moved to
%     "nocleave" section of "allpepdata"
% 1 - glycopeptide only
% 2 - glycopeptide and peptide
% 3 - peptide only
if ismember(input.analyzetgt,[1,2,3])
    glypep = allpepdata.glypep;
    digestable_list = false(size(glypep));
    ptmtyp = allpepdata.ptminfo.uniind;
    % Get PTM markers, e.g. 1 for Met oxidation, 2 for N/O-glycan, 3 for Cys carboxymethylation
    % If a glycopeptide carries a marker "2", then the final library will
    %     contain glycopeptides that carries every member of N-glycan
    %     library
    ptmseq = allpepdata.ptminfo.mod;
    [unityp,ind] = unique(ptmtyp);
    isglycan = false(size(unityp));
    for ii = 1:length(ind)
        isglycan(ii) = any(strfind(ptmseq{ind(ii)},'{'));  % Glycans are recognized by "{"
    end
    glycanind = ind(isglycan);
    switch input.analyzetgt
        case 1  % N-Glycopeptide only
            for ii = 1:length(glypep)
                thisprot = glypep{ii};
                keepglypep = false(size(thisprot));
                % Identify the PTM marker on each glypep, if peptide does
                %     not carry glycan, abandon this one (Case 1: glypep only)
                for jj = 1:length(keepglypep)
                    ptmsites = regexp(thisprot{jj},'{[0-9]+}','match');
                    ptmind = zeros(size(ptmsites));  % Find the markers, check if they are glycans.
                    for k = 1:length(ptmsites)
                        ptmind(k) = str2double(ptmsites{k}(2:end-1));
                    end
                    if any(ismember(ptmind,glycanind))  % PTMs identified on peptide include glycan
                        if sum(ismember(ptmind,glycanind)) == 1
                            keepglypep(jj) = true;
                        end
                    end
                end
                if any(keepglypep)
                    digestable_list(ii) = true;  % If a protein contains no glypep it's rejected
                    glypep{ii} = thisprot(keepglypep);
                end
            end
        case 2 % O-Glycopeptide Only
            
        case 3  % Peptide only
            for ii = 1:length(glypep)
                thisprot = glypep{ii};
                keepglypep = false(size(thisprot));
                for jj = 1:length(keepglypep)
                    ptmsites = regexp(thisprot{jj},'{[0-9]+}','match');
                    ptmind = zeros(size(ptmsites));
                    for k = 1:length(ptmsites)
                        ptmind(k) = str2double(ptmsites{k}(2:end-1));
                    end
                    if ~any(ismember(ptmind,glycanind))  % No glycan PTM found
                        keepglypep(jj) = true;
                    end
                end
                if any(keepglypep)
                    digestable_list(ii) = true;
                    glypep{ii} = thisprot(keepglypep);
                end
            end
    end
    FASTAhead = allpepdata.FASTAhead;
    allpepdata.FASTAhead = FASTAhead(digestable_list);
    allpepdata.glypep = glypep(digestable_list);
    allpepdata.nocleave = [allpepdata.nocleave;FASTAhead(~digestable_list)];
    PTMisglycan = allpepdata.ptminfo.isglycan;
    PTMuniind = allpepdata.ptminfo.uniind;
    PTMglycanind = unique(PTMuniind(PTMisglycan));
    for ii = 1:length(allpepdata.glypep)
        tempprot = allpepdata.glypep{ii};
        keepind = true(size(tempprot));
        for jj = 1:length(tempprot)
            [~,g,~] = breakGlyPep(tempprot{jj});
            ptmindonpep = zeros(size(g));
            for kk = 1:length(g)
                ptmindonpep(kk) = str2double(g(kk).struct(2:end-1));
            end
            if sum(ismember(ptmindonpep,PTMglycanind)) > 1
                keepind(jj) = false;
            end
        end
        allpepdata.glypep{ii} = tempprot(keepind);
    end
    protkeepind = ~cellfun(@isempty,allpepdata.glypep);
    allpepdata.FASTAhead = allpepdata.FASTAhead(protkeepind);
    allpepdata.glypep = allpepdata.glypep(protkeepind);
    allpepdata.nocleave = [allpepdata.nocleave;allpepdata.FASTAhead(~protkeepind)];
end

if ~isempty(statusreporthandles)  % Report progress to GUI, if exists
    statusstr = [statusstr;'Performing PTM fragmentation...'];
    set(statusreporthandles.edit_statusreport, 'String', statusstr);
    pause(0.0001);
end

% Parallel computing significantly saves time, but needs to be set up first
numpw = input.doparacomp;
doparallel = numpw > 1;
activepool = gcp('nocreate');
if doparallel  % is parallel computing toolbox installed?
    scoreoptions.doparallel = true;
    scoreoptions.parallelactive = true;
    if isempty(activepool)
        parpool('local',numpw);
    end
else
    if ~isempty(activepool)
        delete(gcp('nocreate'));
    end
end
if scoreoptions.resumeprevsearch
    % If there is an interrrupted previous search, this function can find
    %     out where it stopped and where to continue.
    todo = analyzeexistingfiles(input,allpepdata);
else
    todo = cell(1,2);
    todo{1} = 1:length(allpepdata.glypep);
    todo{2} = 1:length(allpepdata.glypep);
end
if input.presearch
    if ~isempty(todo{1})
        % Tailored input for prec search
        % Designed for N-glycan only - for now
        % IMPORTANT: Here we assume all N-glycans occupy every site possible
        %     (N-X-S/T, X ~= P)
        NGlycan_pres = {'{n{n{h{h}{h{h}{h}}}}}';'{n{n{h{h}{h{h}}}}}';'{n{n{h{h}{h}}}}';...
            '{n{f}{n{h{h{n}}{h{n}}}}}';'{n{n{h{h{n{h}}}{h}}}}';...
            '{n{f}{n{h{h{n{h}}}{n}{h}}}}';'{n{f}{n{h{h{n{h}}}{h{n{h}}}}}}';'{n{n{h{h{n{h}}}{h{n{h{s}}}}}}}';...
            '{n{n{h{h{n{h{s}}}}{h{n{h{s}}}}}}}';'{n{f}{n{h{h{n{h{s}}}}{h{n{h{s}}}}}}}';...
            '{n{n{h{h{n}}{h{n}{n}}}}}';'{n{n{h{h{n{h{s}}}}{h{n}{n{h}}}}}}';...
            '{n{n{h{h{n{h{s}}}}{h{n}{n{h{s}}}}}}}';'{n{n{h{h{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}';'{n{n{h{h{n}{n}}{h{n}{n}}}}}';...
            '{n{n{h{h{n{h{s}}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}';'{n{n{h{h{n{h{s}}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}';'{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}';...
            '{n{n{h{h{n{h}}}{h{n{h}}}}}}';'{n{f}{n{h{h{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}';'{n{f}{n{h{h{n{h}}}{h{n{h}}{n{h}}}}}}';...
            '{n{n{h{h{n{h}}}{h{n{h}}{n{h}}}}}}'};
        
        % BELOW: for test purposes only, 200+ glycans 
        % NGlycan_pres = {'{n{n{h{h}{h}}}}','{n{n{h{h{n}}{h}}}}','{n{n{h{h{n}}{h{n}}}}}','{n{n{h{h{n{n{s}}}}{h}}}}','{n{n{h{h{n}}{h{n}{n}}}}}','{n{n{h{h{n{n{s}}}}{h{n}}}}}','{n{n{h{h{n}{n}}{h{n}{n}}}}}','{n{n{h{h{n{n}}}{h{n{n{s}}}}}}}','{n{n{h{h{n{n{s}}}}{h{n{n{s}}}}}}}','{n{n{h{h{n}{n}}{n}{h{n}{n}}}}}','{n{n{h{h}{h{h}}}}}','{n{n{h{h{n{h}}}{h}}}}','{n{n{h{h{n{h{s}}}}{h}}}}','{n{n{h{h{n{h}}}{n}{h}}}}','{n{n{h{h{n{h{s}}}}{n}{h}}}}','{n{n{h{h{n{h}}}{h{n}{n}}}}}','{n{n{h{h{n{h{s}}}}{h{n}{n}}}}}','{n{n{h{h{n{n{s}}}}{h{n{h{s}}}}}}}','{n{n{h{h{n{h}}}{n}{h{n}{n}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n}{n}}}}}','{n{n{h{h{n}{n}}{n}{h{n}{n{h}}}}}}','{n{n{h{h{n}{n}}{n}{h{n}{n{h{s}}}}}}}','{n{n{h{h}{h{h}{h}}}}}','{n{n{h{h{n{h}}}{h{h}}}}}','{n{n{h{h{n{h{s}}}}{h{h}}}}}','{n{n{h{h{n{h}}}{h{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{h{n{h{s}}}}}}}','{n{n{h{h{n{h}}}{h{n}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{h{n}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{h{n}{n{h{s}}}}}}}','{n{n{h{h{n{h}}}{n}{h{n}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n}{n{h{s}}}}}}}','{n{n{h{h{n}{n{h}}}{n}{h{n}{n{h}}}}}}','{n{n{h{h{n}{n{h}}}{n}{h{n}{n{h{s}}}}}}}','{n{n{h{h{n}{n{h{s}}}}{n}{h{n}{n{h{s}}}}}}}','{n{n{h{h{h}}{h{h}{h}}}}}','{n{n{h{h{n{h}}}{h{h}{h}}}}}','{n{n{h{h{n{h{s}}}}{h{h}{h}}}}}','{n{n{h{h{n{h{s}}}}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{n}{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n}{n{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n}{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{h}}{h{h}{h{h}}}}}}','{n{n{h{h{n{h}}{n{h}}}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n{h}}{n{h}}}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}{n{h{s}}}}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{n{h{s}}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{n{h{h{n{h}}{n{h}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{n{h{h{n{h}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{n{h{s}}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{n{h{h{h{h}}}{h{h}{h{h}}}}}}','{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}','{n{f}{n{h{h}{h}}}}','{n{f}{n{h{h{n}}{h}}}}','{n{f}{n{h{h{n}}{h{n}}}}}','{n{f}{n{h{h{n{n{s}}}}{h}}}}','{n{f}{n{h{h{n}}{h{n}{n}}}}}','{n{f}{n{h{h{n{n{s}}}}{h{n}}}}}','{n{f}{n{h{h{n}{n}}{h{n}{n}}}}}','{n{f}{n{h{h{n}{n}}{n}{h{n}{n}}}}}','{n{f}{n{h{h}{h{h}}}}}','{n{f}{n{h{h{n{h}}}{h}}}}','{n{f}{n{h{h{n{h{s}}}}{h}}}}','{n{f}{n{h{h{n{h}}}{n}{h}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h}}}}','{n{f}{n{h{h{n{h}}}{h{n}{n}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n}{n}}}}}','{n{n{h{h{n{n{s}}}}{h{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}}{n}{h{n}{n}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n}{n}}}}}','{n{f}{n{h{h{n}{n}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n}{n}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h}{h{h}{h}}}}}','{n{f}{n{h{h{n{h}}}{h{h}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{h}}}}}','{n{f}{n{h{h{n{h}}}{h{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h}}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n}{n{h}}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h{s}}}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}}{h{h}{h}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{h}{h}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n}{n{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n}{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h{s}}}}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n}{n}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n}{n}}}}}','{n{f}{n{h{h{n{n{s}}}}{h{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n}{n}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n}{n}}}}}','{n{f}{n{h{h{n}{n}}{n}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n}}{n}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{h}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{h}}}}}','{n{f}{n{h{h{n{h{s}}}}{h{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h}}}{n}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{h{s}}}}{n}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{h}{h}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{h}{h}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{h}}{n{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{h{s}}}}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{h}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{n}{h{n{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{h{s}}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{f}{h{s}}}}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{n}{h{n{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{n}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}}{n}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n}{n{f}{h{s}}}}{n}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{f}{h{s}}}}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h}}{n{f}{h}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{n}{h{n{f}{h}}{n{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{n}{h{n{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{h{s}}}{n{f}{h{s}}}}{n}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}{n{f}{h{s}}}}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{n}{h{n{f}{h}}{n{f}{h}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h}}}{n}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{n}{h{n{f}{h}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h}}{n{f}{h{s}}}}{n}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{f}{n{h{h{n{f}{h{s}}}{n{f}{h{s}}}}{n}{h{n{f}{h{s}}}{n{f}{h{s}}}}}}}','{n{n{h{h}{n}{h}}}}','{n{n{h{h{n}}{n}{h}}}}','{n{n{h{h{n}}{n}{h{n}}}}}','{n{n{h{h{n}}{n}{h{n}{n}}}}}','{n{n{h{h{n}}{h{h}}}}}','{n{n{h{h{n}}{h{n{h{s}}}}}}}','{n{n{h{h{n}}{h{n{h{s}}}{n}}}}}','{n{n{h{h{n}{n}}{h{n}{n{h}}}}}}','{n{n{h{h{n}{n}}{h{n}{n{h{s}}}}}}}','{n{n{h{h{n}}{h{h}{h}}}}}','{n{n{h{h{n{h}}}{h{n{h{s}}}}}}}'}';
        
        
        % Use the whole peptide library, but replace N-glycans in variable PTM
        %    library with the pre-search specific one.
        allpepdata_pres = allpepdata;
        varptm_pres = allpepdata_pres.varptm;
        ptminfo_pres = allpepdata_pres.ptminfo;
        isNglycan = strcmpi(varptm_pres.aaresidue,'N');
        Nglycanlocator = unique(varptm_pres.varptmuniind(isNglycan));
        ptminfo_keepind = ptminfo_pres.uniind ~= Nglycanlocator;
        ptminfo_pres.uniind = [ptminfo_pres.uniind(ptminfo_keepind);...
            repmat(Nglycanlocator,length(NGlycan_pres),1)];
        ptminfo_pres.mod = [ptminfo_pres.mod(ptminfo_keepind);NGlycan_pres];
        allpepdata_pres.ptminfo = ptminfo_pres;
        % The above steps look for N-glycans, replace the PTM strings and keep
        %     everything else the old way.
        scoreoptions_presearch = scoreoptions;
        scoreoptions_presearch.presearch = 1;
        searchspectra(input,msdata,triggerdata,allpepdata_pres,statusreporthandles,...
            todo{2},input.allfragments_presearch,scoreoptions_presearch);
        
        % Collect and analyze data
        [keepwhichprotein,~] = analyzepresearchdata_test(input,triggerdata,allpepdata_pres,...
            scoreoptions_presearch,5);
        % Using the results, we now know which proteins have a higher chance of
        %     being detected. Build a much smaller search space - use the
        %     original variable PTM and selected proteins.
        allpepdata_filtered = allpepdata;
        allpepdata_filtered.FASTAhead = allpepdata.FASTAhead(keepwhichprotein);
        allpepdata_filtered.glypep = allpepdata.glypep(keepwhichprotein);
        scoreoptions.presearch = 2;
        [~,f,~] = fileparts(input.outputfname);
        % Save the result of pre-search for the record.
        save(fullfile(input.outputdir,[f,'_intercept.mat']),'allpepdata','keepwhichprotein','allpepdata_pres','allpepdata_filtered');
        save(fullfile(input.outputdir,[f,'_filteredglypep.mat']),'allpepdata_filtered');
        
        pause(3);
        if isempty(allpepdata_filtered.FASTAhead) || isempty(allpepdata_filtered.glypep)
            error('No glycoprotein found after Pre-Search');
        end
        searchspectra(input,msdata,triggerdata,allpepdata_filtered,statusreporthandles,...
            todo{2},input.allfragments,scoreoptions);
        fbioptions.nummindiagionfragmode = 2;
        fbioptions.allpepdata = allpepdata_filtered;
        findbestisomers(input,msdata,input.allfragments,scoreoptions,fbioptions);
    end
else
    if ~isempty(todo{2})
        searchspectra(input,msdata,triggerdata,allpepdata,statusreporthandles,...
            todo{2},scoreoptions);
    end
end
end