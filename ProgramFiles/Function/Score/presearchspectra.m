function presearchresults = presearchspectra(input,msdata,triggerdata,allpepdata,statusreporthandles,...
    groupstodo,scoreoptions)
% PRESEARCHSPECTRA: match MS2 spectra against user provided glycopeptide
%     candidates.
% 
% Syntax:
% searchspectra(input,msdata,triggerdata,numpw,allpepdata,...
%     statusreporthandles,scoreoptions)
% 
% Input:
% input: structure. Information necessary for scoring. Same as parent
%     function's.
% msdata: structure. Experiment data. Same as parent function's.
% triggerdata: structure. Scan number association in HCD trigger
%     experiment. Fields are:
%         SASSO: n x m numerical array. Scan numbers and parent MS1 scan
%             retention time. For each row, the first n - 2 scans'
%             fragmentation methods are defined in field "colnames", the
%             second to last number is the parent MS1 scan number, the last
%             number is the retention time of it.
%         colnames: 1 x n cell array of strings. The fragmentation methods
%             of the MS2 scans.
% allpepdata: structure. The information in digested glycopeptide file. See
%     DIGESTFILEANALY for details.
% statusreporthandles: structure. The handle of scoregui. This is the
%     channel to the textbox which provides realtime status update.
% scoreoptions: structure. Settings and parameters for scoring. The fields
%     are :
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     analyzefragmode            n x 1 cell array of              The fragmentation modes to be analyzed.
%                                                      strings       
%     assemblegpmode           Double                                 Controls working mode of function ASSEMBLEGP, See the
%                                                                                                      function's document for detail.
%     doparallel                          Logical                                  Enable parallel computing mode or not.
%     userpepiontyp                 n x 1 cell array of              Peptide ion types to be considered for each fragmentation mode.
%                                                      strings.
%     (Below are identical to the fields of "input" in parent function, see SCOREALLSPECTRA for detail.)
%     ms1tol
%     ms1tolunit
%     ms2tol
%     ms2tolunit
%     minmaxmz
%     maxlag
%     cutoffmed
%     fracmax
%     fragnum
%     selectpeak
%     numprotperworker
%     monosacislabile
%     simultaneousfrag
%     addoxoniumion
%     isHCDtrigger
%     maxstublen
%     presearch
%     presearchtopprotnum
%     resumeprevsearch
% ---------------------------------------------------------------------------------------------------------------------------------------------------
% 
% Output:
% Scoring results are saved as .mat files. Each file contains 2 fields:
%     "scoreintdata" and "result". "scoreintdata" is a structure containing
%     user settings and intermediate data. These data are useful when
%     browsing the results. "result" is a 1 x n structure where scoring
%     results were stored. See SAVESCORERESULTS for details about
%     "scoreintdata". See MATCH1BY1 for details about "result".
% 
% Note:
% N/A
%
% Example:
% N/A. Use provided sample data to perform a complete scoring process, set
%     breaking points to see how this function works.
%
% Children function: 
% N/A
% 
% See Also:
% SCOREALLSPECTRA THEOPTMFRAG FINDMS1MATCH MATCHSPECTRA SAVESCORERESULTS

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

if ~isempty(statusreporthandles)
    statusstr = get(statusreporthandles.edit_statusreport, 'String');
end

% Experiment data
scannum = msdata.scannum;
charge = msdata.charge;
precmz = msdata.precursormz;
precmass   = (precmz - 1.007825032).*charge;
doparallel = scoreoptions.doparallel;

% Result storage
outputdir = input.outputdir;
outputfilename = input.outputfname;
% ptmstruct = allpepdata.ptminfo.mod;
% ptmfragopt.monosacislabile = scoreoptions.monosacislabile;

% Generate all theoretical fragments of PTMs. Later they will be combined
%     with peptide fragments to build glycopeptide fragments.
% if ~doparallel
%     [PTM.ptmseq,PTM.ptmmass,...
%         PTM.ptmmassfix,PTM.ptmtype,...
%         PTM.ptmfragsto,PTM.ptmfragstublen] = ...
%         theoptmfrag(ptmstruct,scoreoptions.fragnum(:,2),scoreoptions.analyzefragmode,ptmfragopt);
% else  % SPMD
%     if scoreoptions.doparallel && ~scoreoptions.parallelactive
%         parpool('local',input.doparacomp);
%     end
%     % Note: MATLAB has a timeout feature on parallel computing environment.
%     %     Therefore the existence of the environment should be checked
%     %     regularly.
%     ptmstruct_dist = distributed(ptmstruct);
%     spmd
%         [ptmseq_dist,ptmmass_dist,ptmmassfix_dist,ptmtype_dist,ptmfragsto_dist,ptmfragstublen_dist] = ...
%             theoptmfrag(getLocalPart(ptmstruct_dist),scoreoptions.fragnum(:,2),...
%             scoreoptions.analyzefragmode,ptmfragopt);
%     end
%     % For definition of fields please check the document of THEOPTMFRAG
%     PTM.ptmseq = [ptmseq_dist{:}];
%     PTM.ptmmass = [ptmmass_dist{:}];
%     PTM.ptmmassfix = [ptmmassfix_dist{:}];
%     PTM.ptmtype = [ptmtype_dist{:}];
%     PTM.ptmfragsto = [ptmfragsto_dist{:}];
%     PTM.ptmfragstublen = [ptmfragstublen_dist{:}];
%     clear *_dist;
% end

%% Step 4: Start data analysis
% To maximize the efficiency and minimize resource consumption, jobs are
% distributed to each worker in candidate - spectrum pairs, rather than
% protein - msdata. The redundant information distributed to workers is
% now eliminated.

% Find out the spectra that have a proper precursor ion mass first, then
% match against a glycopeptide candidate.

% Identify the spectra to be searched.
if scoreoptions.isHCDtrigger
    % HCD trigger is special because the same ion may be fragmented
    % multiple times. In this case only one search is needed.
    SASSO = triggerdata.SASSO;
    SASSO(SASSO(:,2) == 0,:) = [];  % HCD THAT TRIGGERED NOTHING IS IGNORED!
    triggerdata.SASSO = SASSO;
    ms1searchscannum = SASSO(:,1);
else
    ms1searchscannum = scannum(msdata.mslvl == 2 &...
        ismember(msdata.fragmode,scoreoptions.analyzefragmode));
end

% Proteins are analyzed in batches, the number of proteins to be processed
% each time is decided by user setting.
numprotperbatch = ceil(input.doparacomp * scoreoptions.numprotperworker);
if numprotperbatch < 1
    errordlg({'Workload allocation error.',...
        'Check setting "Each worker analyze __ proteins each time"'},...
        'Check input');
    return
end
totalprotnum = length(allpepdata.glypep);
for i = 1:ceil(totalprotnum/numprotperbatch)
    if ismember(i,groupstodo)
        protindse = [(i-1)*numprotperbatch+1,min(i*numprotperbatch,length(allpepdata.glypep))];
        protbatch = allpepdata.glypep(protindse(1):protindse(2));
        FASTAbatch = allpepdata.FASTAhead(protindse(1):protindse(2));
        protind = protindse(1):protindse(2);
        if ~doparallel
            [matched_gpcomp,matched_gpmass,matched_scanind] = ...
                findms1match(precmass,scannum,protbatch,protind,input.allfragments_presearch,...
                ms1searchscannum,scoreoptions);
        else
            dist_prot = distributed(protbatch);
            dist_protind = distributed(protind);
            spmd
                [dist_matched_gpcomp,dist_matched_gpmass,dist_matched_scanind] = ...
                    findms1match(precmass,scannum,getLocalPart(dist_prot),...
                    getLocalPart(dist_protind),input.allfragments_presearch,ms1searchscannum,scoreoptions);
            end
                matched_gpcomp = [dist_matched_gpcomp{:}];
                matched_gpmass = [dist_matched_gpmass{:}];
                matched_scanind = [dist_matched_scanind{:}];
            clear dist_*
        end
        if ~isempty(statusreporthandles)
            analytime = toc;
            statusstr = [['PRESEARCH: Group ',num2str(i),': Protein No.',num2str(protindse(1)),...
                ' to ',num2str(protindse(2)),' starteded, time = ',num2str(analytime),' sec.'];...
                statusstr];
            set(statusreporthandles.edit_statusreport, 'String', statusstr);
            pause(0.0001);
        end
        result = presearch_matchspectrA(msdata,matched_gpcomp,matched_gpmass,...
            matched_scanind,protbatch,FASTAbatch,triggerdata,protind,input.allfragments_presearch,scoreoptions);
        pause(5);
        savescoreresults(result,scoreoptions,input,input.allfragments_presearch,allpepdata,triggerdata,...
            outputdir,outputfilename,protindse,'presearchsave');
        if ~isempty(statusreporthandles)
            analytime = toc;
            statusstr = [['PRESEARCH: Group ',num2str(i),': Protein No.',num2str(protindse(1)),...
                ' to ',num2str(protindse(2)),' finished, time = ',num2str(analytime),' sec.'];...
                [num2str(totalprotnum-protindse(2)),' remaining.'];'Saving results...';statusstr];
            set(statusreporthandles.edit_statusreport, 'String', statusstr);
            pause(0.0001);
        end
    end
end
end