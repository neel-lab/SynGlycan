function [SASSO,colnames] = buildm2sasso(input,options)
% BUILDM2SASSO: build ms2 supplemental activatation association
%
% Syntax:
%   [SASSO,colnames] = buildm2sasso(input,options)
%
% Input:
%
% Output:
%   SASSO structure: [T , F1 , F2 , ... , MS1_T, MS1_RT]
%     T: trigger scan number
%     F1, F2, ...: follow up scans triggered by T
%     MS1_T: parent MS1 scan of T
%     MS1_RT: retenteion time of MS1_T
%   colnames: Fragmentation modes that correspond to all F's in SASSO
%
% Examples:
%   [SASSO,colnames] = buildm2sasso(input,options)
%
%See also: PREPROCESSGUI
%
% NOTE: This program assumes that 1 trigger scan is associated with only 1
%       set of follow up scans.
%

workingmode = options.mode;
trigger = options.triggermode;
follow = options.analymode;
if ischar(follow)
    follow = {follow};
end
colnames = [trigger,reshape(follow,1,[])];
switch upper(workingmode)
    case 'PID'
        SASSO = [];
        % one fragmentation triggers multiple subsequent scans
        msdata = input;
        scannum = msdata.scannum;
        retime = msdata.retime;
        precscannum = cell2mat(msdata.precursorScanNum);
        fragmode = msdata.fragmode;
        mslevel = msdata.mslvl;
        triggersn = scannum(strcmpi(fragmode,trigger));
        triggerprecsn = precscannum(strcmpi(fragmode,trigger));
        for ii = 1:length(follow)
            if ~any(strcmpi(fragmode,follow{ii})) || strcmpi(follow{ii},trigger)
                follow{ii} = [];
            end
        end
        follow = follow(~cellfun('isempty',follow));
        colnames = [trigger,reshape(follow,1,[])];
        followsn = cell(size(follow));
        followprecsn = cell(size(follow));
        ms1sn = scannum(mslevel == 1);
        precsnstat = 1;
        for i = 1:length(follow)
            followsn{i} = scannum(strcmpi(fragmode,follow{i}));
            tempfollowprecscans = precscannum(strcmpi(fragmode,follow{i}));
            % DESISION: is precursor scan number right?
            if any(ismember(tempfollowprecscans,ms1sn))
                precsnstat = 2;  % triggered scan has ms1 as prec. scan, need to find the HCD trigger scan
            end
            followprecsn{i} = tempfollowprecscans;
        end
        if precsnstat == 1  % triggered scan's precursor scan numbers are normal
            [~,longerind] = max(cellfun(@length,followsn));
            starter = followsn{longerind};
            starteref = followprecsn{longerind};
            otherfragseq = setdiff(1:length(followsn),longerind);
            followers = [triggersn,followsn(otherfragseq)];
            followerefs = [triggersn,followprecsn(otherfragseq)];
            columnseq = [2,1,otherfragseq+1];
            isreusable = false(size(followers));
            isreusable(1) = true;  % trigger is reusable
            [subSASSO,followers_left,followerefs_left] = getreftable(starter,starteref,followers,...
                followerefs,columnseq,isreusable);
            SASSO = [SASSO;subSASSO];
            whichreusable = find(isreusable);
            followers_left{whichreusable} = setdiff(followers_left{whichreusable},...
                subSASSO(:,whichreusable));
            followerefs_left{whichreusable} = followers_left{whichreusable};
            leftovers_num = cellfun(@length,followers_left);  % exclude HCDs
            leftovers_num(whichreusable) = 0;
            while any(leftovers_num > 0)
                [~,choosestarter] = max(leftovers_num);
                newstarter = followers_left{choosestarter};
                newstarteref = followerefs_left{choosestarter};
                newfollowers_left = [0,followers_left(setdiff(1:length(followers_left),choosestarter))];
                newfollowerefs_left = [0,followerefs_left(setdiff(1:length(followerefs_left),choosestarter))];
                newcolumnseq = [columnseq(choosestarter+1),...
                    columnseq(1),setdiff(columnseq,columnseq([1,choosestarter+1]))];
                % [new base column, prev. exhausted base col, remaining col, ...]
                % new base column will never be trigger
                isreusable = [false,isreusable(setdiff(1:length(followers_left),choosestarter))];
                [subSASSO,followers_left,followerefs_left] = getreftable(newstarter,newstarteref,...
                    newfollowers_left,newfollowerefs_left,newcolumnseq,isreusable);
                SASSO = [SASSO;subSASSO];
                whichreusable = find(isreusable);
                followers_left{whichreusable} = setdiff(followers_left{whichreusable},...
                    subSASSO(:,newcolumnseq(whichreusable+1)));
                followerefs_left{whichreusable} = followers_left{whichreusable};
                leftovers_num = zeros(size(followers_left));
                for i = 1:length(followers_left)
                    leftovers_num(i) = sum(followers_left{i} > 0);
                end
                leftovers_num(whichreusable) = 0;  % exclude HCDs
            end
            triggersn_notmatched = setdiff(triggersn,SASSO(:,1));
            SASSOaddon = zeros(length(triggersn_notmatched),size(SASSO,2));
            SASSOaddon(:,1) = triggersn_notmatched;
            SASSO = [SASSO;SASSOaddon];
            %  at the end of SASSO there's a HCD only region, in case it's useful
        elseif  precsnstat == 2  % need to rebuild precsnstat
            if ms1sn(end) == scannum(end)
                ms1section = [ms1sn(1:end-1),ms1sn(2:end)-1];
            else
                ms1section = [ms1sn(1:end),[ms1sn(2:end)-1;scannum(end)]];
            end
            ms1section = ms1section(diff(ms1section,1,2) > 0,:);  % delete section that holds no MS2 scans
            precmz = msdata.precursormz;
            triggerprecmz = precmz(strcmpi(fragmode,trigger));
            for i = 1:length(followprecsn)
                thesefollowsn = followsn{i};
                thesefollowprecsn = followprecsn{i};
                thesefollowprecmz = precmz(strcmpi(fragmode,follow{i}));
                for j = 1:length(thesefollowprecsn)
                    thisfollowsn = thesefollowsn(j);
                    inms1section = ms1section((thisfollowsn > ms1section(:,1)) &...
                        (thisfollowsn <= ms1section(:,2)),:);
                    insectiontriggersnind = (triggersn > inms1section(1)) &...
                        (triggersn <= inms1section(2));
                    insectiontriggermz = triggerprecmz(insectiontriggersnind);
                    insectiontriggerscannum = triggersn(insectiontriggersnind);
                    identicalprecmzind = abs(insectiontriggermz - thesefollowprecmz(j)) < 0.0001;
                    if any(identicalprecmzind)
                        fixedprecsn = insectiontriggerscannum(identicalprecmzind);
                        thesefollowprecsn(j) = fixedprecsn(1);
                    end
                end
                followprecsn{i} = thesefollowprecsn;
            end
            [~,longerind] = max(cellfun(@length,followsn));
            starter = followsn{longerind};
            starteref = followprecsn{longerind};
            otherfragseq = setdiff(1:length(followsn),longerind);
            followers = [triggersn,followsn(otherfragseq)];
            followerefs = [triggersn,followprecsn(otherfragseq)];
            columnseq = [2,1,otherfragseq+1];
            isreusable = false(size(followers));
            isreusable(1) = true;  % trigger is reusable
            [subSASSO,followers_left,followerefs_left] = getreftable(starter,starteref,followers,...
                followerefs,columnseq,isreusable);
            SASSO = [SASSO;subSASSO];
            whichreusable = find(isreusable);
            followers_left{whichreusable} = setdiff(followers_left{whichreusable},...
                subSASSO(:,whichreusable));
            followerefs_left{whichreusable} = followers_left{whichreusable};
            leftovers_num = cellfun(@length,followers_left);  % exclude HCDs
            leftovers_num(whichreusable) = 0;
            while any(leftovers_num > 0)
                [~,choosestarter] = max(leftovers_num);
                newstarter = followers_left{choosestarter};
                newstarteref = followerefs_left{choosestarter};
                newfollowers_left = [0,followers_left(setdiff(1:length(followers_left),choosestarter))];
                newfollowerefs_left = [0,followerefs_left(setdiff(1:length(followerefs_left),choosestarter))];
                newcolumnseq = [columnseq(choosestarter+1),...
                    columnseq(1),setdiff(columnseq,columnseq([1,choosestarter+1]))];
                % [new base column, prev. exhausted base col, remaining col, ...]
                % new base column will never be trigger
                isreusable = [false,isreusable(setdiff(1:length(followers_left),choosestarter))];
                [subSASSO,followers_left,followerefs_left] = getreftable(newstarter,newstarteref,...
                    newfollowers_left,newfollowerefs_left,newcolumnseq,isreusable);
                SASSO = [SASSO;subSASSO];
                whichreusable = find(isreusable);
                followers_left{whichreusable} = setdiff(followers_left{whichreusable},...
                    subSASSO(:,newcolumnseq(whichreusable+1)));
                followerefs_left{whichreusable} = followers_left{whichreusable};
                leftovers_num = zeros(size(followers_left));
                for i = 1:length(followers_left)
                    leftovers_num(i) = sum(followers_left{i} > 0);
                end
                leftovers_num(whichreusable) = 0;  % exclude HCDs
            end
            triggersn_notmatched = setdiff(triggersn,SASSO(:,1));
            SASSOaddon = zeros(length(triggersn_notmatched),size(SASSO,2));
            SASSOaddon(:,1) = triggersn_notmatched;
            SASSO = [SASSO;SASSOaddon];
        end
        SASSO = sortrows(SASSO,1);
        ms1_ms1retime = zeros(size(SASSO,1),2);
        for i = 1:size(SASSO,1)
            thisms1 = triggerprecsn(triggersn == SASSO(i,1));
            thisms1rt = retime(scannum == thisms1);
            ms1_ms1retime(i,:) = [thisms1,thisms1rt];
        end
        SASSO = [SASSO,ms1_ms1retime];
    case 'MULTIFILEPID'
        
end
end












