function [cleanerseq,specialoptions] = getglydressing(GLYseq)
% GETGLYDRESSING: analyze glycan sequence (glycan + customization),
%     return separated information
%
% Syntax:
% [cleanerseq,specialoptions] = getglydressing(GLYseq)
%
% Input:
% GLYseq: string, IUPAC sequence of glycan, sequence may contain
%     customization information: anomer, adduct, curly bracket, etc.
%
% Output:
% cleanerseq: string. IUPAC sequence of glycan, deprived of customization.
% specialoptions: structure. Customization information of glycan. Fieldnames
%     are the name of options, corresponding values are n x 2 cell array,
%     1st column contains value of the option, 2nd is the serial number of
%     monosac. it applies to. In most cases, each cell contains only 1
%     monosac., except for "CURLYBRACKET".
%
% Note:
% This function handles glycan only.
% This function is specifically designed to handle ambiguous structures
%     represented by being enclosed by a pair of curly brackets. All other
%     options are handled through GETSEQOPTIONS.
% Ambiguous structures (curly bracket items) must be placed at the
%     beginning of the main structure.
% 
%
% Example:
% [cleanerseq,specialoptions] = getglydressing('{Neu5Ac,Neu5Ac}{Man[Man]Man[GlcNAc[GlcNAc]Man]Man(-adduct "Na^{+}")GlcNAc()GlcNAc()UDP(-char)}')
%
% cleanerseq =
%
% Man[Man]Man[GlcNAc[GlcNAc]Man]Man()GlcNAc()GlcNAc()UDP()
%
%
% specialoptions =
%
%     CURLYBRACKET: {'Neu5Ac,Neu5Ac'  [1 2 3 4 5 6 7 8 9 10]}
%           ADDUCT: {'Na^{+}'  [7]}
%             CHAR: {[]  [10]}
%
% Children function:
% N/A
% 
% See also:
% GETSEQOPTIONS, SPLIT
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
%


%% Locate content inside quotation mark & block these info
[optcontent,optcontentstart,optcontentend] = regexp(GLYseq,...
    DrawGlycanPara.regexp_optionvalue,'match','start','end');
tempGPseq = GLYseq;
if ~isempty(optcontent)  % block customization info ("noise")
    for i = 1:length(optcontentstart)
        tempGPseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i) - ...
            optcontentstart(i) + 1);
    end
end
specialoptions = [];
%%  Curly bracket handling
indtemp = strfind(tempGPseq,'{');
% After replacing local option values with whichspaces  there should be no
% curly bracket anywhere except those related to ambiguous structures
if any(indtemp)
    % Below is the standard curlybracket counter used everywhere
    levelindex = zeros(2,length(tempGPseq));
    indtemp2 = strfind(tempGPseq,'}');
    if length(indtemp) ~= length(indtemp2)
        errordlg('Unsymmetrical curly brackets.','Check input - Ambiguous structure');
        return
    elseif length(indtemp) ~= 2 || length(indtemp2) ~= 2
        errordlg('Only 2 pairs of curly brackets allowed.','Check input - Ambiguous structure');
        return
    end
    
    levelindex(1,indtemp) = 1;
    levelindex(1,indtemp2) = -1;
    levelindex(2,1) = 1;
    for i = 2:size(levelindex,2)
        levelindex(2,i) = levelindex(2,i-1) + levelindex(1,i);
    end
    numsection = length(indtemp)/2+1;
    sectionpos = zeros(numsection,2);
    sectionend = find(levelindex(2,:) == 0)';
    sectionpos(:,2) = [sectionend(1:numsection-1);sectionend(end)];
    sectionpos(1,1) = 1;
    sectionpos(2:end,1) = sectionpos(1:end-1,2)+1;
    cbitems = cell(size(sectionpos,1)-1,2);
    for i = 1:size(cbitems,1)
        cbitems{i,1} = GLYseq(sectionpos(i,1)+1:sectionpos(i,2)-1);
    end
    
    mainstructGPseq = tempGPseq(sectionpos(end,1):sectionpos(end,2));
    % This part is complicated - originally designed for allowing bracket
    % to cover only part of the main structure. This function was
    % cancelled, but it can be used to handle minor input errors.
    indtempmain = strfind(mainstructGPseq,'{');
    levelindex = zeros(1,length(mainstructGPseq));
    indtemp2main = strfind(mainstructGPseq,'}');
    levelindex(indtempmain) = 1;
    levelindex(indtemp2main) = -1;
    opencbrep = [];
    cbpair = [];
    ind = 1;
    while ind <= length(levelindex)
        if levelindex(ind) == 1
            opencbrep = [opencbrep;ind];
        elseif levelindex(ind) == -1
            cbpair = [cbpair;[opencbrep(end),ind]];
            opencbrep(end) = [];
        end
        ind = ind + 1;
    end
    cbpair = sortrows(cbpair,1);
    cbpos = regexp(mainstructGPseq,'[{}]');
    mainstructGPseq = regexprep(mainstructGPseq,'[{}]','');
    if strcmp(mainstructGPseq(end),' ')
        mainstructGPseq(end) = '';
    end
    [~,MonoIndex,~,~] = Split(mainstructGPseq);
    for i = 1:length(MonoIndex)
        MonoIndex(i) = MonoIndex(i) + sum(MonoIndex(i) > cbpos);
    end
    cbcontent = cell(size(cbpair,1),1);
    for i = 1:size(cbpair,1)
        cbcontent{i} = find(MonoIndex < cbpair(i,2) & MonoIndex >= cbpair(i,1));
    end
    for i = 1:size(cbitems,1)
        cbitems{i,2} = cbcontent{i};  % curly bracket coverage area, shown using monosac. index number
    end
    
    specialoptions.CURLYBRACKET = cbitems;
    cleanGLYseq = GLYseq(sectionpos(end,1):sectionpos(end,2));
    cleanGLYseq([indtempmain,indtemp2main]) = '';
else
    cleanGLYseq = GLYseq;
end
%%  End of curly bracket handling

%% other contents
[otheroptions,cleanerseq] = getseqoptions(cleanGLYseq,'G');
% Combine local options in sequence with directly injected options
if ~isempty(otheroptions)
    otheroptionnames = fieldnames(otheroptions);
    for i = 1:length(otheroptionnames)
        if isfield(specialoptions,otheroptionnames{i})
            specialoptions.(otheroptionnames{i}) = [specialoptions.(otheroptionnames{i});...
                otheroptions.(otheroptionnames{i})];
        else
            specialoptions.(otheroptionnames{i}) = otheroptions.(otheroptionnames{i});
        end
    end
end
end
