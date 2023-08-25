function [SASSO,colnames] = buildm2sasso(msdata,options)
% BUILDM2SASSO: build MS2 supplemental activatation association between
%     scans
%
% Syntax:
% [SASSO,colnames] = buildm2sasso(msdata,options)
%
% Input:
% msdata: structure. .mat MS data file.
% options: structure. Field names are:
%     mode: string. Available options are:
%         'PID': product ion dependent
%         'MULTIPLEPID': (unavailable now) a place holder.
%     triggermode: string. The fragmentation method that triggered
%         supplemental scans.
%     analymode: 1 x n cell array of strings. Fragmentation methods that
%         are triggered.
%
% Output:
% SASSO: m x n numerical array. m equals to the number of trigger scans. n
%     equals to number of "analymode" + 3
%     SASSO structure in each row: [T , F1 , F2 , ... , MS1_T, MS1_RT]
%         T: trigger scan number
%         F1, F2, ...: supplemental scans triggered by T
%         MS1_T: parent MS1 scan of T
%         MS1_RT: retenteion time of MS1_T
% colnames: 1 x n cell array of strings. Fragmentation modes of
%     trigger and supplemental scans.
%
% NOTE: This program assumes that 1 trigger scan is associated with only 1
%     set of follow up scans. In the event of 1 scan triggers multiple
%     sets, for example scan a triggers (b1, c1) and (b2,c2), where b and c
%     are different fragmentation methods, SASSO will be shown as:
%     [a, b1, c1]
%     [a, b2, c2]
%
% Examples:
%   [SASSO,colnames] = buildm2sasso(msdata,options)
%
% Children function:
% N/A
%
% See Also:
% PREPROCESSGUI

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

workingmode = options.mode;
trigger = options.triggermode;
follow = options.analymode;
if iscell(trigger)
    trigger = trigger{1};
end
if ischar(follow)
    follow = {follow};
end
colnames = [trigger,reshape(follow,1,[])];
switch upper(workingmode)
    case 'PID'
        % One fragmentation triggers multiple subsequent scans
        scannum = msdata.scannum;
        retime = msdata.retime;
        precscannum = cell2mat(msdata.precursorScanNum);
        fragmode = msdata.fragmode;
        mslevel = msdata.mslvl;
        triggersn = scannum(strcmpi(fragmode,trigger));
        triggerprecsn = precscannum(strcmpi(fragmode,trigger));
        for ii = 1:length(follow)  % follow: supplemental (triggered) scans
            if ~any(strcmpi(fragmode,follow{ii})) || strcmpi(follow{ii},trigger)
                follow{ii} = [];  % Only consider fragmentation modes that is present
            end
        end
        follow = follow(~cellfun('isempty',follow));
        colnames = [trigger,reshape(follow,1,[])];
        followsn = cell(size(follow));
        followprecsn = cell(size(follow));
        ms1sn = scannum(mslevel == 1);
        precsnstat = 1;
        for ii = 1:length(follow)
            followsn{ii} = scannum(strcmpi(fragmode,follow{ii}));
            tempfollowprecscans = precscannum(strcmpi(fragmode,follow{ii}));
            % DESISION: is precursor scan number right?
            % Method: In processed .mat file, the field "precursorScanNum"
            %     contains parent scan number for each MS2 scan. For
            %     triggered scans, their parent scan should be the triggering
            %     scan, whose MS level = 2.
            if precsnstat == 1 && any(ismember(tempfollowprecscans,ms1sn))
                precsnstat = 2;
            end
            followprecsn{ii} = tempfollowprecscans;
        end
        followsnnums = cellfun(@length,followsn);
        if precsnstat == 1  % precursor scan numbers are correct
            SASSO = [];
            for ii = 1:length(triggersn)
                temptriggersn = triggersn(ii);
                tempSASSO = [];
                followscannum = cell(size(follow));
                for jj = 1:length(follow)
                    tempfollowscannum = followsn{jj}(followprecsn{jj} == temptriggersn);
                    if ~isempty(tempfollowscannum)
                        followscannum{jj} = tempfollowscannum;
                    end
                end
                if ~any(~cellfun(@isempty,followscannum))
                    
                elseif length(unique(cellfun(@length,followscannum))) == 1
                    for jj = 1:length(follow)
                        tempfollowsn = followscannum{jj};
                        if jj == 1
                            tempSASSO = repmat(temptriggersn,length(tempfollowsn),length(follow) + 1);
                        end
                        tempSASSO(:,jj + 1) = tempfollowsn(:);
                    end
                else
                    errordlg('Trigger scan unequal for different frag. mode.','Unable to Proceed')
                end
                SASSO = [SASSO;tempSASSO];
            end
            % At the end of SASSO there's a HCD only region, these scans
            %     are kept for potential uses
            triggernf = setdiff(triggersn,SASSO(:,1));  % Scans of trigger fragmentation mode but triggered nothing
            SASSO(followsnnums(1)+1:followsnnums(1) + length(triggernf),1) = triggernf;
        elseif  precsnstat == 2  % Triggered scans have wrong parent scan, fix this issue
            SASSO = zeros(0,length(follow) + 1);
            % Assuming that between 2 MS1 scans, the parent scan of all MS2
            %     scans within is the beginning MS1 scan. These 2 MS1 scans
            %     form a "section".
            % To rebuild triggered scan's parent scan, within each section
            %     first find the triggering scans, then for each of them
            %     find triggered scans with identical precursor m/z. Since
            %     triggering scans is always in front of triggered scans,
            %     the scans are naturally sorted.
            if ms1sn(end) == scannum(end)
                ms1section = [ms1sn(1:end-1),ms1sn(2:end)-1];
            else
                ms1section = [ms1sn(1:end),[ms1sn(2:end)-1;scannum(end)]];
            end
            precmz = msdata.precursormz;
            triggerprecmz = precmz(strcmpi(fragmode,trigger));
            followprecmz = cell(size(follow));
            for ii = 1:length(follow)
                followprecmz{ii} = precmz(strcmpi(fragmode,follow{ii}));
            end
            ms1section = ms1section(diff(ms1section,1,2) > 0,:);  % delete section that holds no MS2 scans
            for ii = 1:size(ms1section,1)
                trigger_irind = triggersn > ms1section(ii,1) & triggersn <= ms1section(ii,2);  % ir = in range
                trigger_ir = triggersn(trigger_irind);
                trigger_irpmz = triggerprecmz(trigger_irind);
                follow_ir = followsn{1}(followsn{1} > ms1section(ii,1) & followsn{1} <= ms1section(ii,2));
                follow_irpmz = followprecmz{1}(followsn{1} > ms1section(ii,1) & followsn{1} <= ms1section(ii,2));
                tempSASSO = [];
                if any(follow_irpmz == 0) && any(follow_irpmz > 0)  % mix of 0s and > 0s
                    errordlg('Precursor Ion Mass Has Zero and non-Zeros','Unable to Proceed');
                    return
                elseif any(follow_irpmz == 0) && ~any(follow_irpmz > 0)  % all 0
                    trigger_irdiff = find(diff(trigger_ir) > 1);
                    trigger_irdiff = trigger_irdiff(:);
                    followse = [];
                    for jj = 1:length(trigger_irdiff)
                        followse = [followse;[trigger_ir(trigger_irdiff(jj)) + 1,trigger_ir(trigger_irdiff(jj) + 1) - 1]];
                    end
                    followse = [followse;[trigger_ir(end) + 1,ms1section(ii,2)]];
                    triggerse = trigger_ir([[1;trigger_irdiff + 1],[trigger_irdiff;length(trigger_ir)]]);
                    for jj = 1:size(followse,1)
                        currenttrigger = triggerse(jj,1):triggerse(jj,2);
                        numtrigger = length(currenttrigger);
                        numfollow = followse(jj,2) - followse(jj,1) + 1;
                        numfollowpertrigger = numfollow / numtrigger / length(follow);
                        temptempSASSO = zeros(numfollow/length(follow),length(follow) + 1);
                        if mod(numfollowpertrigger,1) == 0  % integer meaning each trigger has equal number of follower
                            for kk = 1:length(follow)
                                follow_ir = followsn{kk}(followsn{kk} >= followse(jj,1) & followsn{kk} <= followse(jj,2));
                                if kk == 1
                                    for ll = 1:numtrigger
                                        temptempSASSO((ll - 1)*numfollowpertrigger + 1:ll*numfollowpertrigger,1) = ...
                                            repmat(currenttrigger(ll),numfollowpertrigger,1);
                                    end
                                end
                                temptempSASSO(:,kk + 1) = follow_ir;
                            end
                        end
                        tempSASSO = [tempSASSO;temptempSASSO];
                    end
                else  % all > 0
                    for jj = 1:length(trigger_irpmz)
                        temptempSASSO = [];
                        for kk = 1:length(follow)
                            follow_ir = followsn{kk}(followsn{kk} > ms1section(ii,1) & ...
                                followsn{kk} <= ms1section(ii,2));
                            follow_irpmz = followprecmz{kk}(followsn{kk} > ms1section(ii,1) & ...
                                followsn{kk} <= ms1section(ii,2));
                            isprecmzmatch = find(follow_irpmz == trigger_irpmz(jj));
                            if any(isprecmzmatch)
                                thisfollowmatch = follow_ir(isprecmzmatch);
                                temptempSASSO(:,kk + 1) = thisfollowmatch(:);
                                if kk == 1
                                    temptempSASSO(:,1) = repmat(trigger_ir(jj),length(thisfollowmatch),1);
                                end
                            end
                        end
                        tempSASSO = [tempSASSO;temptempSASSO];
                    end
                end
                SASSO = [SASSO;tempSASSO];
            end
            triggernf = setdiff(triggersn,SASSO(:,1));  % Scans of trigger fragmentation mode but triggered nothing
            SASSO(followsnnums(1)+1:followsnnums(1)+length(triggernf),1) = triggernf;
        end
        
        % Sort by triggering scans first, then separate triggered and
        %     untriggered, put untriggered at the end.
        SASSO = sortrows(SASSO,1);
        triggeronly = zeros(size(SASSO,1),1);
        for ii = 1:length(follow)
            triggeronly = triggeronly | SASSO(:,ii+1) == 0;
        end
        triggeronly = logical(triggeronly);
        SASSO = [SASSO(~triggeronly,:);SASSO(triggeronly,:)];
        for ii = 1:size(SASSO,1)
            thisms1 = triggerprecsn(triggersn == SASSO(ii,1));
            thisms1rt = retime(scannum == thisms1);
            SASSO(ii,length(colnames)+1:length(colnames)+2) = [thisms1,thisms1rt];
        end
        
    case 'MULTIFILEPID'
        % INTENTIONALLY LEFT BLANK
end
end