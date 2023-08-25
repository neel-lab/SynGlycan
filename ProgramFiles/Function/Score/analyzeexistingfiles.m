function todo = analyzeexistingfiles(input,allpepdata)
% ANALYZEEXISTINGFILES: analyze existing scoring files, then return the
%     jobs to be done after an interrupted previous search.
%
% Syntax:
% todo = analyzeexistingfiles(input)
%
% Input:
% input: structure. The umberlla varaible holding all input information
%     including scoring options and experiment data. The detail of its
%     fields can be found in function SCOREALLSPECTRA.
%
% Output:
%  todo: 1 x 2 cell array of n x 1 numerical array. The groups yet to be
%      analyzed. 1st element is for pre-search (if enabled), 2nd element is
%      for main search. The members of each group is calculated by
%      input parameters, therefore it is fixed throughout the entire
%      analysis. Downstream functions will use these group serial numbers
%      to determin which group to analyze.
%
% Note:
% It is imperative to maintain the consistency of scoring parameters,
%     including storage location of result files, MS1 and MS2 tolerence
%     settings, denoising powers, etc. Make sure these parameters remains
%     the same as the previous search.
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% N/A

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

totalprotnum = length(allpepdata.glypep);
numprotperbatch = ceil(input.doparacomp * input.numprotperworker);
numgroups = ceil(totalprotnum/numprotperbatch);
outputfilename = input.outputfname;
[~,f,e] = fileparts(outputfilename);
filenames = cell(numgroups,4);  % {presearchname,mainsearchname}
for ii = 1:numgroups
    protindse = [(ii-1)*numprotperbatch+1,min(ii*numprotperbatch,totalprotnum)];
    filenames{ii,1} = [f,'_presearch','_',num2str(protindse(1)),'_',num2str(protindse(2)),e];
    filenames{ii,2} = [f,'_',num2str(protindse(1)),'_',num2str(protindse(2)),e];
    % If file too large, data is divided and stored in different files.
    % If everything goes right
    % special case of 1: this is a regexp pattern
    filenames{ii,3} = [f,'_presearch','_',num2str(protindse(1)),'_',...
        num2str(protindse(2)),'_Part_[0-9]+of[0-9]+',e];
    % special case of 2: this is a regexp pattern too
    filenames{ii,4} = [f,'_',num2str(protindse(1)),'_',...
        num2str(protindse(2)),'_Part_[0-9]+of[0-9]+',e];
end
todo_presearch = [];
todo_mainsearch = [];
tgtfolderdir = dir(input.outputdir);
existingdatafilenames = {tgtfolderdir.name};
existingdatafilenames = existingdatafilenames(~[tgtfolderdir.isdir]);
if input.presearch == 1
    presearchfilebeencreated = false(size(filenames,1),1);
    for ii = 1:size(filenames,1)
        case1 = ismember(filenames{ii,1},existingdatafilenames);
        % when result is divided
        case3 = false;  % Initialization
        allcase3 = ~cellfun(@isempty,regexp(existingdatafilenames,filenames{ii,3}));
        if any(allcase3)
            case3filenames = existingdatafilenames(allcase3);
            [~,firstcase3filename,~] = fileparts(case3filenames{1});
            groupsizestringstart = strfind(firstcase3filename,'of');
            groupsizestringstart = groupsizestringstart(end);
            groupsize = str2double(firstcase3filename(groupsizestringstart+2:end));
            if sum(allcase3) == groupsize
                case3 = true;
            end
        end
        presearchfilebeencreated(ii) = case1 || case3;
    end
    if any(~presearchfilebeencreated)
        todo_presearch = find(~presearchfilebeencreated);
    end
    mainsearchfilebeencreated = false(size(filenames,1),1);
    for ii = 1:size(filenames,1)
        case2 = ismember(filenames{ii,2},existingdatafilenames);
        % when result is divided
        case4 = false;  % Initialization
        allcase4 = ~cellfun(@isempty,regexp(existingdatafilenames,filenames{ii,4}));
        if any(allcase4)
            case4filenames = existingdatafilenames(allcase4);
            [~,firstcase4filename,~] = fileparts(case4filenames{1});
            groupsizestringstart = strfind(firstcase4filename,'of');
            groupsizestringstart = groupsizestringstart(end);
            groupsize = str2double(firstcase4filename(groupsizestringstart+2:end));
            if sum(allcase4) == groupsize
                case4 = true;
            end
        end
        mainsearchfilebeencreated(ii) = case2 || case4;
    end
    if any(~mainsearchfilebeencreated)
        todo_mainsearch = find(~mainsearchfilebeencreated);
    end
elseif input.presearch == 0  % Direct search
    mainsearchfilebeencreated = false(size(filenames,1),1);
    for ii = 1:size(filenames,1)
        case2 = ismember(filenames{ii,2},existingdatafilenames);
        % when result is divided
        case4 = false;  % Initialization
        allcase4 = ~cellfun(@isempty,regexp(existingdatafilenames,filenames{ii,4}));
        if any(allcase4)
            case4filenames = existingdatafilenames(allcase4);
            [~,firstcase4filename,~] = fileparts(case4filenames{1});
            groupsizestringstart = strfind(firstcase4filename,'of');
            groupsizestringstart = groupsizestringstart(end);
            groupsize = str2double(firstcase4filename(groupsizestringstart+2:end));
            if sum(allcase4) == groupsize
                case4 = true;
            end
        end
        mainsearchfilebeencreated(ii) = case2 || case4;
    end
    if any(~mainsearchfilebeencreated)
        todo_mainsearch = find(~mainsearchfilebeencreated);
    end
end
todo = {todo_presearch,todo_mainsearch};
end