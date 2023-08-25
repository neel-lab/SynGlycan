function specialoptions = getdgfrag(PGM,theofrag,matchresult)
% GETDGFRAG: generate fragment marker text and position for DrawGlycan-SNFG
% to use
%
% Syntax:
% specialoptions = getdgfrag(PGM,theofrag,matchresult)
%
% Input:
% PGM: [PGM{1},PGM{2},PGM{3}] = breakGlyPep(SGP)
% theofrag: Theoretical fragments of candidate (glyco)peptide.
% matchresult: result returned by CALCITHSCORE, must contain field
%     "peakmatchindex".
%
% Output:
% specialoptions: n x 1 cell array, each element contains fragmentation
%     marker info for each PTM, cells are empty for non-glycan PTMs or no
%     fragmentation occured. Also peptide fragmentation markers are placed in
%     the first cell. The whole variable is designed to be readable by
%     DrawGlycan-SNFG.
%
% NOTE:
% "glybreak" is glycosidic bond break only, peptide is intact,
%     because these fragments has a special iontype format. "pepbreak"
%     contains peptide break only, even if there is co-fragmentation, the
%     glycan break is not included
%
% Example:
% N/A
%
% Children function:
% FINDGLYBONDBREAK
%


[glybondbreakpos,unifrag] = findglybondbreak(PGM,theofrag,matchresult);
p = PGM{1};
g = PGM{2};
m = PGM{3};
allptmpos = [];
specialoptions = {};
if ~isempty(g)
    allptmpos = [allptmpos,[g.pos;zeros(size(g))]];
end
if ~isempty(m)
    allptmpos = [allptmpos,[m.pos;ones(size(m))]];
end  % ptmpos & identity
if ~isempty(allptmpos)
    [~,ind] = sort(allptmpos(1,:));
    allptmpos = allptmpos(:,ind);
    specialoptions = cell(size(allptmpos,2),1);
    for ii = 1:size(allptmpos,2)
        if allptmpos(2,ii) == 1
            specialoptions{ii}.CHAR = {[],1};
        end
    end
end
pepbreakind = [];
pepbreakinfo = {};
glybreakind = [];
glybreakinfo = {};
if ~isempty(p)
    aanum = length(p.pos);
else
    aanum = 0;
end
if ~isempty(unifrag)
    for ii = 1:length(unifrag)
        unitindex = unifrag(ii).unitindex;
        thispepfragunitind = unitindex{1};
        thisglyfragunitind = unitindex{2};
        thistyp = -1;
        % thistyp:
        % default = -1
        % Gly B-ion = 0
        % Gly Y-ion = 1
        % Pep abc-ion (w/ glycan co-fragmentation) = 2
        % pep xyz-ion (w/ glycan co-fragmentation) = 3
        % Pep abc-ion (w/o glycan co-fragmentation) = 4
        % pep xyz-ion (w/o glycan co-fragmentation) = 5
        iontype = unifrag(ii).type;
        if ismember(iontype(1),{'a','b','c'})  % pep/co-frag
            thistyp = 2;
            if ~any(strfind(iontype,'-'))
                thistyp = 4;
            end
        elseif ismember(iontype(1),{'x','y','z'})  % pep/co-frag
            thistyp = 3;
            if ~any(strfind(iontype,'-'))
                thistyp = 5;
            end
        elseif ismember(iontype(1),{'A','B','C'})  % N/A, but keep it here
        elseif ismember(iontype(1),{'X','Y','Z'})  % N/A, but keep it here
        elseif strcmpi(iontype(1),'-') && ismember(iontype(2),{'a','b','c'})  % "-b" means gly frag
            thistyp = 0;
        elseif strcmpi(iontype(1),'-') && ismember(iontype(2),{'x','y','z'})  % "-y" means gly frag
            thistyp = 1;
        elseif ismember(iontype(2),{'A','B','C'})  % N/A, but keep it here
        elseif ismember(iontype(2),{'X','Y','Z'})  % N/A, but keep it here
        end
        if strcmpi(iontype(1),'-')  % B/Y
            if thistyp == 0  % a/b/c
                ionstring = ['B',num2str(size(unifrag(ii).unitindex{2},1))];
            elseif thistyp == 1
                ionstring = ['Y',num2str(size(unifrag(ii).unitindex{2},1))];
            end
        else  % co-frag/"none"
            iontypecont = strsplit(iontype,'-');
            ionstring = iontypecont{1};
        end
        if ~isempty(thispepfragunitind) && length(thispepfragunitind) < aanum  % pep frag
            if ismember(thistyp,[2,4]) % a/b/c - cofrag/no cofrag
                pepbreakind = [pepbreakind;[max(thispepfragunitind),thistyp]];
                if thistyp == 2
                    pepbreakinfo = [pepbreakinfo;[ionstring,'*']];
                else
                    pepbreakinfo = [pepbreakinfo;ionstring];
                end
            elseif ismember(thistyp,[3,5])  % x/y/z - cofrag/no cofrag
                pepyfrag = min(thispepfragunitind);
                pepbreakind = [pepbreakind;[pepyfrag,thistyp]];
                if thistyp == 3
                    pepbreakinfo = [pepbreakinfo;[ionstring,'*']];
                else
                    pepbreakinfo = [pepbreakinfo;ionstring];
                end
            end
        end
        if ~isempty(thisglyfragunitind)
            if thistyp == 0 % A/B/C
                glybreakmarkpos = min(thisglyfragunitind(:,1));
                if ~isempty(glybreakmarkpos)
                    glybreakind = [glybreakind;[glybreakmarkpos,thistyp,ii]];
                    glybreakinfo = [glybreakinfo;ionstring];
                end
            elseif thistyp == 1  % X/Y/Z
                glybreakmarkpos = glybondbreakpos{ii,1};
                if ~isempty(glybreakmarkpos)
                    for jj = 1:length(glybreakmarkpos)
                        glybreakind = [glybreakind;[glybreakmarkpos(jj),thistyp,ii]];
                        glybreakinfo = [glybreakinfo;ionstring];
                    end
                else
                    glybreakind = [glybreakind;[1,thistyp,ii]];
                    glybreakinfo = [glybreakinfo;ionstring];
                end
            end
        else
            if thistyp == 1  % Y0
                glybreakmarkpos = glybondbreakpos{ii,1};
                if isempty(glybreakmarkpos)
                    glybreakind = [glybreakind;[1,thistyp,ii]];
                    glybreakinfo = [glybreakinfo;ionstring];
                end
            end
        end
    end
    %% KEEP ANNOTATION TO MINIMUM
    % frag. marker may overlay, avoid overlays, when unavoidable,
    % combine them to a new marker
    if ~isempty(pepbreakind)
        % PEPTIDE - STAGE 1 - NOT DISTINGUISHING CO-FRAG
        oripepbreakind = pepbreakind;
        oripepbreakinfo = pepbreakinfo;
        pepbreakind = [];
        pepbreakinfo = {};
        [~,~,unitind] = unique(oripepbreakind(:,1));
        for ii = 1:max(unitind)
            temppepbreakind = oripepbreakind(unitind == ii,:);
            temppepbreakinfo = oripepbreakinfo(unitind == ii);
            % abc-type pep frag
            tempuniinfostring = '';
            pepbreaktyp2ind = temppepbreakind(:,2) == 2;
            if any(pepbreaktyp2ind)
                tempwritepepbreakind = temppepbreakind(pepbreaktyp2ind,:);
                unipepbreaktyp2info = unique(temppepbreakinfo(pepbreaktyp2ind));
                for jj = 1:length(unipepbreaktyp2info)
                    tempuniinfostring = [tempuniinfostring,unipepbreaktyp2info{jj}(1),','];
                end
                tempuniinfostring = [tempuniinfostring(1:end-1),unipepbreaktyp2info{1}(2:end),'; '];
            end
            pepbreaktyp4ind = temppepbreakind(:,2) == 4;
            if any(pepbreaktyp4ind)
                tempwritepepbreakind = temppepbreakind(pepbreaktyp4ind,:);
                unipepbreaktyp4info = unique(temppepbreakinfo(pepbreaktyp4ind));
                for jj = 1:length(unipepbreaktyp4info)
                    tempuniinfostring = [tempuniinfostring,unipepbreaktyp4info{jj}(1),','];
                end
                tempuniinfostring = [tempuniinfostring(1:end-1),unipepbreaktyp4info{1}(2:end)];
            else
                tempuniinfostring = tempuniinfostring(1:end-2);
            end
            if ~isempty(tempuniinfostring)
                pepbreakind = [pepbreakind;[tempwritepepbreakind(1,1),2]];
                pepbreakinfo = [pepbreakinfo;tempuniinfostring];
            end
            % xyz-type pep frag
            tempuniinfostring = '';
            pepbreaktyp3ind = temppepbreakind(:,2) == 3;
            if any(pepbreaktyp3ind)
                tempwritepepbreakind = temppepbreakind(pepbreaktyp3ind,:);
                unipepbreaktyp3info = unique(temppepbreakinfo(pepbreaktyp3ind));
                for jj = 1:length(unipepbreaktyp3info)
                    tempuniinfostring = [tempuniinfostring,unipepbreaktyp3info{jj}(1),','];
                end
                tempuniinfostring = [tempuniinfostring(1:end-1),unipepbreaktyp3info{1}(2:end),'; '];
            end
            pepbreaktyp5ind = temppepbreakind(:,2) == 5;
            if any(pepbreaktyp5ind)
                tempwritepepbreakind = temppepbreakind(pepbreaktyp5ind,:);
                unipepbreaktyp5info = unique(temppepbreakinfo(pepbreaktyp5ind));
                for jj = 1:length(unipepbreaktyp5info)
                    tempuniinfostring = [tempuniinfostring,unipepbreaktyp5info{jj}(1),','];
                end
                tempuniinfostring = [tempuniinfostring(1:end-1),unipepbreaktyp5info{1}(2:end)];
            else
                tempuniinfostring = tempuniinfostring(1:end-2);
            end
            if ~isempty(tempuniinfostring)
                pepbreakind = [pepbreakind;[tempwritepepbreakind(1,1),3]];
                pepbreakinfo = [pepbreakinfo;tempuniinfostring];
            end
        end
    end
    % GLYCAN
    if ~isempty(glybreakind)
        keepwhich = true(size(glybreakind,1),1);
        ii = 1;
        occupation = zeros(max(glybreakind(:,1)),1);
        while ii <= size(glybreakind,1)
            if keepwhich(ii)
                identicalind = find(glybreakind(:,3) == glybreakind(ii,3));
                prevmarkerpos = glybreakind(1:ii,1);
                [~,tempind] = setdiff(glybreakind(identicalind,1),prevmarkerpos);
                keepwhich(identicalind) = false;
                if any(tempind)  % empty marker place available, use the first empty one
                    placeind = identicalind(tempind(1));
                else  % must use a prev occupied place, just use the first one
                    [~,tempind2] = min(occupation(glybreakind(identicalind,1)));
                    placeind = identicalind(tempind2);  % this part find the least occupied
                    % space to place the marker
                end
                keepwhich(placeind) = true;
                occupation(glybreakind(placeind,1)) = occupation(glybreakind(placeind,1)) + 1;
            end
            ii = identicalind(end) + 1;
        end
        % The above lines is to make sure each frag ion gets one marker, but they
        % still can be combined further, based on marker position.
        tempgbind = glybreakind(keepwhich,:);
        tempgbinfo = glybreakinfo(keepwhich);
        glybreakind = [];
        glybreakinfo = {};
        [~,~,ind] = unique(tempgbind(:,1:2),'rows');
        for ii = 1:max(ind)
            tempind = find(ind == ii);
            glybreakind = [glybreakind;tempgbind(tempind(1),1:2)];
            tempbreakinfo = unique(tempgbinfo(tempind));
            tempbreakinfostring = tempbreakinfo{1}(1);
            for jj = 1:length(tempbreakinfo)
                tempbreakinfostring = [tempbreakinfostring,tempbreakinfo{jj}(2:end),','];
            end
            glybreakinfo = [glybreakinfo;tempbreakinfostring(1:end-1)];
        end
    end
    %% BUILD DRAWGLYCAN COMPATIBLE INPUT
    pepbreak = [];
    if ~isempty(pepbreakind)
        if any(pepbreakind(:,2) == 2)  % abc
            temp = pepbreakind(pepbreakind(:,2) == 2,:);
            tempinfo = pepbreakinfo(pepbreakind(:,2) == 2,:);
            temppepbreakcell = cell(size(temp,1),2);
            for ii = 1:size(temp,1)
                temppepbreakcell{ii,1} = tempinfo{ii};
                temppepbreakcell{ii,2} = temp(ii,1) ;
            end
            pepbreak.N = temppepbreakcell;
        end
        if any(pepbreakind(:,2) == 3)  % xyz
            temp = pepbreakind(pepbreakind(:,2) == 3,:);
            tempinfo = pepbreakinfo(pepbreakind(:,2) == 3,:);
            temppepbreakcell = cell(size(temp,1),2);
            for ii = 1:size(temp,1)
                temppepbreakcell{ii,1} = tempinfo{ii};
                temppepbreakcell{ii,2} = temp(ii,1) - 1;
            end
            pepbreak.C = temppepbreakcell;
        end
    end
    if ~isempty(allptmpos)
        glybreak = cell(size(allptmpos,2),1);
        glycanlengths = 0;
        if ~isempty(g)
            glycanlengths = [g.len];
        end
        msnumag = glycanlengths*triu(ones(length(glycanlengths)) - eye(length(glycanlengths)),0);
        ptmposshift = find(allptmpos(2,:) == 0);
        % because in DrawGlycan PTM is treated as glycan, frag marker might be
        % placed in a spot reserved for non-glycan PTM
        if ~isempty(glybreakind)
            for ii = 1:size(glybreakind,1)
                tempglybreakind = glybreakind(ii,1);
                whichgly = find(tempglybreakind > msnumag,1,'last');  % find which glycan it belongs to
                if isempty(whichgly)  % that's the first glycan
                    whichgly = 0;  % this number is an offset value to correctify frag info
                    fragmarkeroffset = 0;
                else
                    fragmarkeroffset = msnumag(whichgly);
                end
                if glybreakind(ii,2) == 1  % Y-ion
                    glybreak{ptmposshift(whichgly)} = [glybreak{ptmposshift(whichgly)};...
                        {glybreakind(ii,1) - fragmarkeroffset,'R',glybreakinfo{ii}}];
                elseif glybreakind(ii,2) == 0  % B-ion
                    glybreak{ptmposshift(whichgly)} = [glybreak{ptmposshift(whichgly)};...
                        {glybreakind(ii,1) - fragmarkeroffset,'NR',glybreakinfo{ii}}];
                end
            end
        end
        for ii = 1:length(specialoptions)
            if ~isempty(glybreak{ii})
                specialoptions{ii}.NR = {};
                specialoptions{ii}.R = {};
            end
        end
        for ii = 1:length(glybreak)
            if ~isempty(glybreak{ii})
                isNR = strcmpi(glybreak{ii}(:,2),'NR');
                tempNRfrag = [glybreak{ii}(isNR,3),glybreak{ii}(isNR,1)];
                specialoptions{ii}.NR = [specialoptions{ii}.NR;tempNRfrag];
                isR = strcmpi(glybreak{ii}(:,2),'R');
                tempRfrag = [glybreak{ii}(isR,3),glybreak{ii}(isR,1)];
                specialoptions{ii}.R = [specialoptions{ii}.R;tempRfrag];
            end
        end
    end
    if ~isempty(pepbreak)
        if isfield(pepbreak,'N')
            specialoptions{1}.pepN = pepbreak.N;
        end
        if isfield(pepbreak,'C')
            specialoptions{1}.pepC = pepbreak.C;
        end
    end
end
end