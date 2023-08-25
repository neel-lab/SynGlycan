function [displaydata,displaydataind] = rebuilddisptable(scoredata,displaydataind,...
    displayoptions,usercustom,sevenfilterbys)
% Isomers are deleted by sets (deleting 1 isomer means deleting it from all
%     its fragmentation modes.
filterbyfragmode = sevenfilterbys{1};
% which method: combies/every hit/each frag individually,...
filterbyglypeponly = sevenfilterbys{2};  % TRUE/FALSE
filterbyes = sevenfilterbys{3};  % TRUE/FALSE
filterbyescutoff = sevenfilterbys{4};  % 0 ~ 1
filterbyspecqual = sevenfilterbys{5};  % TRUE/FALSE
filterbyprecion = sevenfilterbys{6};  % TRUE/FALSE
filterbyuserdec = sevenfilterbys{7};  % TRUE/FALSE
scoreintdata = displayoptions.scoreintdata;
combineisomer = displayoptions.combineisomer;
isHCDtrigger = scoreintdata.scoreoptions.isHCDtrigger;
quality3 = usercustom.Quality3;
allenscore = [scoredata.Enscore];
allsgpseq = {scoredata.SGP};
allfragmode = {scoredata.Fragmode};

numfragmethods = 1;
if isHCDtrigger
    allfragmethod = scoreintdata.colnames;
    numfragmethods = length(allfragmethod);
else
    allfragmethod = scoreintdata.scoreoptions.analyzefragmode;
end


if combineisomer
    displaydataind_height = sum(cellfun(@length,displaydataind(:,1)));
    expanddisplaydataind = zeros(displaydataind_height,numfragmethods);
    cellindex = zeros(displaydataind_height,1);
    currentind = 1;
    for ii = 1:size(displaydataind,1)
        for jj = 1:numfragmethods
            tempdisplaydataind = displaydataind{ii,jj};
            numnewline = size(tempdisplaydataind,1);
            expanddisplaydataind(currentind:currentind + numnewline - 1,jj) = ...
                displaydataind{ii,jj};
        end
        cellindex(currentind:currentind + numnewline - 1) = ii;
        currentind = currentind + numnewline;
    end
else
    expanddisplaydataind = displaydataind;
end

if strcmpi(filterbyfragmode,'CombiES')  % must be HCD trigger
    [CombiES,~] = calculateCombiES(scoredata,expanddisplaydataind,...
        allfragmethod);
    keeprow = false(size(expanddisplaydataind,1),1);
    for ii = 1:size(expanddisplaydataind,1)
        tempdisplaydataind = expanddisplaydataind(ii,:);
        isomerCombiES = CombiES(ii);
        representind = tempdisplaydataind(tempdisplaydataind > 0);
        representind = representind(1);
        sgpseq = allsgpseq{representind};
        keepthis = checkquality(sgpseq,isomerCombiES,quality3,representind,...
            filterbyes,filterbyescutoff,filterbyglypeponly,filterbyspecqual,...
            filterbyprecion,filterbyuserdec);
        if ~keepthis
            continue
        end
        keeprow(ii) = true;
    end
elseif strcmpi(filterbyfragmode,'Everyone')
    for ii = 1:size(expanddisplaydataind,1)
        for jj = 1:size(expanddisplaydataind,2)
            displaydataind_backup = expanddisplaydataind(ii,jj);
            expanddisplaydataind(ii,jj) = 0;
            if displaydataind_backup > 0
                keepthis = checkquality(allsgpseq{displaydataind_backup},...
                    allenscore(displaydataind_backup),...
                    quality3,displaydataind_backup,...
                    filterbyes,filterbyescutoff,filterbyglypeponly,filterbyspecqual,...
                    filterbyprecion,filterbyuserdec);
                if ~keepthis
                    continue
                end
                expanddisplaydataind(ii,jj) = displaydataind_backup;
            end
        end
    end
    keeprow = ~all(expanddisplaydataind == 0,2);
elseif strcmpi(filterbyfragmode,'Each frag. individually')
    fragmodes = scoreoptions.analyzefragmode;  % This is the sequence of the multiple ES input.
    if scoreintdata.scoreoptions.isHCDtrigger  % Trigger experiment
        [~,ind1] = sort(allfragmethod);
        [~,ind2] = sort(fragmodes);
        [~,ind3] = sort(ind1);
        ind = ind2(ind3);
        filterbyescutoff = filterbyescutoff(ind);
        for ii = 1:size(expanddisplaydataind,1)
            for jj = 1:size(expanddisplaydataind,2)
                displaydataind_backup = expanddisplaydataind(ii,jj);
                expanddisplaydataind(ii,jj) = 0;
                if displaydataind_backup > 0
                    keepthis = checkquality(allsgpseq{displaydataind_backup},...
                        allenscore(displaydataind_backup),...
                        quality3,displaydataind_backup,...
                        filterbyes,filterbyescutoff(jj),filterbyglypeponly,filterbyspecqual,...
                        filterbyprecion,filterbyuserdec);
                    if ~keepthis
                        continue
                    end
                    expanddisplaydataind(ii,jj) = displaydataind_backup;
                end
            end
        end
    else
        for ii = 1:size(expanddisplaydataind,1)
            displaydataind_backup = expanddisplaydataind(ii);
            expanddisplaydataind(ii) = 0;
            if displaydataind_backup > 0
                tempfragmode = allfragmode{displaydataind_backup};
                tempfilterbyescutoff = filterbyescutoff(ismember(fragmodes,tempfragmode));
                keepthis = checkquality(allsgpseq{displaydataind_backup},...
                    allenscore(displaydataind_backup),...
                    quality3,displaydataind_backup,...
                    filterbyes,tempfilterbyescutoff,filterbyglypeponly,filterbyspecqual,...
                    filterbyprecion,filterbyuserdec);
                if ~keepthis
                    continue
                end
                expanddisplaydataind(ii) = displaydataind_backup;
            end
        end
    end
    keeprow = ~all(expanddisplaydataind == 0,2);
else
    % filter specified frag mode only
    if scoreintdata.scoreoptions.isHCDtrigger
        methodind = find(ismember(allfragmethod,filterbyfragmode));
        if isempty(methodind)
            errordlg('Fragmentation method not found.','Fragmentation method not found.');
        else
            for ii = 1:size(expanddisplaydataind,1)
                displaydataind_backup = expanddisplaydataind(ii,methodind);
                expanddisplaydataind(ii,methodind) = 0;
                if displaydataind_backup > 0
                    keepthis = checkquality(allsgpseq{displaydataind_backup},...
                        allenscore(displaydataind_backup),...
                        quality3,displaydataind_backup,...
                        filterbyes,filterbyescutoff,filterbyglypeponly,filterbyspecqual,...
                        filterbyprecion,filterbyuserdec);
                    if ~keepthis
                        continue
                    end
                    expanddisplaydataind(ii,methodind) = displaydataind_backup;
                end
            end
        end
    else
        for ii = 1:size(expanddisplaydataind,1)
            displaydataind_backup = expanddisplaydataind(ii);
            if displaydataind_backup > 0
                tempfragmode = allfragmode{displaydataind_backup};
               if strcmpi(tempfragmode,filterbyfragmode)
                   keepthis = checkquality(allsgpseq{displaydataind_backup},...
                       allenscore(displaydataind_backup),...
                       quality3,displaydataind_backup,...
                       filterbyes,filterbyescutoff,filterbyglypeponly,filterbyspecqual,...
                       filterbyprecion,filterbyuserdec);
                   if ~keepthis
                       expanddisplaydataind(ii) = 0;
                   end
               end
            end
        end
    end
    keeprow = ~all(expanddisplaydataind == 0,2);
end    
if combineisomer
    newexpanddisplaydataind = expanddisplaydataind(keeprow,:);
    newcellindex = cellindex(keeprow);
    [~,~,ind] = unique(newcellindex);
    displaydataind = cell(max(ind),numfragmethods);
    for ii = 1:max(ind)
        for jj = 1:numfragmethods
            tempdisplaydataind = newexpanddisplaydataind(ind == ii,jj);
            displaydataind{ii,jj} = tempdisplaydataind;
        end
    end
else
    displaydataind = expanddisplaydataind(keeprow,:);
end
displaydata = {};
if ~isempty(displaydataind)    
    [displaydata,~,~,~] = builddisptable(scoredata,displaydataind,displayoptions,usercustom);
end
end

function keepthis = checkquality(sgpseq,enscore,quality3,representind,...
    filterbyes,filterbyescutoff,filterbyglypepyn,filterbyspecqual,...
    filterbyprecion,filterbyuserdec)
keepthis = false;
if filterbyes
    if enscore < filterbyescutoff
        return
    end
end
if filterbyglypepyn
    if ~any(strfind(sgpseq,'{')) || ~any(strfind(sgpseq,'}'))
        return
    end
end
if filterbyspecqual && ~quality3(representind,1)
    return
end
if filterbyprecion && ~quality3(representind,2)
    return
end
if filterbyuserdec && ~quality3(representind,3)
    return
end
keepthis = true;
end


% function [displaydata,displaydataind] = rebuilddisptable(scoredata,displaydataind,...
%     displayoptions,usercustom,sevenfilterbys)
% % Isomers are deleted by sets (deleting 1 isomer means deleting it from all
% %     its fragmentation modes.
% filterbyfragmode = sevenfilterbys{1};
% % which method: combies/every hit/each frag individually,...
% filterbyglypeponly = sevenfilterbys{2};  % TRUE/FALSE
% filterbyes = sevenfilterbys{3};  % TRUE/FALSE
% filterbyescutoff = sevenfilterbys{4};  % 0 ~ 1
% filterbyspecqual = sevenfilterbys{5};  % TRUE/FALSE
% filterbyprecion = sevenfilterbys{6};  % TRUE/FALSE
% filterbyuserdec = sevenfilterbys{7};  % TRUE/FALSE
% combineisomer = displayoptions.combineisomer;
% scoreintdata = displayoptions.scoreintdata;
% quality3 = usercustom.Quality3;
% allenscore = [scoredata.Enscore];
% allsgpseq = {scoredata.SGP};
% allfragmode = {scoredata.Fragmode};
% 
% if strcmpi(filterbyfragmode,'CombiES')  % must be HCD trigger
%     SASSO_method = scoreintdata.colnames;
%     [CombiES,~] = calculateCombiES(scoredata,displaydataind,...
%         SASSO_method);
%     numfragmethods = length(SASSO_method);
%     if combineisomer  % CombiES - HCD trigger - combine isomer
%         for i = 1:size(displaydataind,1)  % each isomer group
%             for j = 1:length(displaydataind{i,1})  % each isomer
%                 displaydataind_backup = zeros(1,numfragmethods);
%                 for k = 1:numfragmethods
%                     displaydataind_backup(k) = displaydataind{i,k}(j);
%                     displaydataind{i,k}(j) = 0;
%                 end
%                 representind = displaydataind_backup(displaydataind_backup > 0);
%                 representind = representind(1);
%                 isomerCombiES = CombiES{i}(j);
%                 if filterbyes && isomerCombiES < filterbyescutoff
%                     continue
%                 end
%                 sgpseq = allsgpseq(representind);
%                 keepthis = checkquality(sgpseq,quality3,representind,...
%                     filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                 if ~keepthis
%                     continue
%                 end
%                 for k = 1:numfragmethods
%                     displaydataind{i,k}(j) = displaydataind_backup(k);
%                 end
%             end
%         end
%     else  % CombiES - HCD trigger - don't combine isomer
%         for i = 1:size(displaydataind,1)  % each isomer
%             displaydataind_backup = displaydataind(i,:);
%             displaydataind(i,:) = zeros(1,numfragmethods);
%             representind = displaydataind_backup(displaydataind_backup > 0);
%             representind = representind(1);
%             isomerCombiES = CombiES(i);
%             if filterbyes && isomerCombiES < filterbyescutoff
%                 continue
%             end
%             sgpseq = allsgpseq(representind);
%             keepthis = checkquality(sgpseq,quality3,representind,...
%                 filterbyglypepyn,filterbyspecqual,filterbyprecion,filterbyuserdec);
%             if ~keepthis
%                 continue
%             end
%             displaydataind(i,:) = displaydataind_backup;
%         end
%     end
% elseif strcmpi(filterbyfragmode,'Every hit')
%     if scoreintdata.scoreoptions.isHCDtrigger
%         SASSO_method = scoreintdata.colnames;
%         numfragmethods = length(SASSO_method);
%         if combineisomer  % Every hit - HCD trigger - combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer group
%                 for j = 1:length(displaydataind{i,1})  % each isomer
%                     for k = 1:numfragmethods
%                         displaydataind_backup = displaydataind{i,k}(j);
%                         displaydataind{i,k}(j) = 0;
%                         if displaydataind_backup > 0
%                             if filterbyes
%                                 tempenscore = allenscore(displaydataind_backup);
%                                 if tempenscore < filterbyescutoff
%                                     continue
%                                 end
%                             end
%                             sgpseq = allsgpseq(displaydataind_backup);
%                             keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                                 filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                             if ~keepthis
%                                 continue
%                             end
%                             displaydataind{i,k}(j) = displaydataind_backup;
%                         end
%                     end
%                 end
%             end
%         else  % Every hit - HCD trigger - don't combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer group
%                 for j = 1:numfragmethods
%                     displaydataind_backup = displaydataind(i,j);
%                     displaydataind(i,j) = 0;
%                     if displaydataind_backup > 0
%                         if filterbyes
%                             tempenscore = allenscore(displaydataind_backup);
%                             if tempenscore < filterbyescutoff
%                                 continue
%                             end
%                         end
%                         sgpseq = allsgpseq(displaydataind_backup);
%                         keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                             filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                         if ~keepthis
%                             continue
%                         end
%                         displaydataind(i,j) = displaydataind_backup;
%                     end
%                 end
%             end
%         end
%     else
%         if combineisomer  % Every hit - regular - combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer
%                 for j = 1:length(displaydataind{i})
%                     displaydataind_backup = displaydataind{i}(j);
%                     displaydataind{i}(j) = 0;
%                     if filterbyes
%                         tempenscore = allenscore(displaydataind_backup);
%                         if tempenscore < filterbyescutoff
%                             continue
%                         end
%                     end
%                     sgpseq = allsgpseq(displaydataind_backup);
%                     keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                         filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                     if ~keepthis
%                         continue
%                     end
%                     displaydataind{i}(j) = displaydataind_backup;
%                 end
%             end
%         else  % Every hit - regular - don't combine isomer
%             for i = 1:size(displaydataind,1)
%                 displaydataind_backup = displaydataind(i);
%                 displaydataind(i) = 0;
%                 if filterbyes
%                     tempenscore = allenscore(displaydataind_backup);
%                     if tempenscore < filterbyescutoff
%                         continue
%                     end
%                 end
%                 sgpseq = allsgpseq(displaydataind_backup);
%                 keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                     filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                 if ~keepthis
%                     continue
%                 end
%                 displaydataind(i) = displaydataind_backup;
%             end
%         end
%     end
% elseif strcmpi(filterbyfragmode,'Each frag. individually')
%     fragmodes = scoreoptions.analyzefragmode;  % This is the sequence of the multiple ES input.
%     if scoreintdata.scoreoptions.isHCDtrigger  % Trigger experiment
%         SASSO_method = scoreintdata.colnames;
%         numfragmethods = length(SASSO_method);
%         if combineisomer  % Each frag. individually - HCD trigger - combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer group
%                 for j = 1:length(displaydataind{i,1})  % each isomer
%                     for k = 1:numfragmethods
%                         displaydataind_backup = displaydataind{i,k}(j);
%                         displaydataind{i,k}(j) = 0;
%                         if displaydataind_backup > 0
%                             if filterbyes
%                                 tempenscore = allenscore(displaydataind_backup);
%                                 thisfragmethod = allfragmode(displaydataind_backup);
%                                 tempescutoff = filterbyescutoff(ismember(fragmodes,thisfragmethod));
%                                 if tempenscore < tempescutoff
%                                     continue
%                                 end
%                             end
%                             sgpseq = allsgpseq(displaydataind_backup);
%                             keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                                 filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                             if ~keepthis
%                                 continue
%                             end
%                             displaydataind{i,k}(j) = displaydataind_backup;
%                         end
%                     end
%                 end
%             end
%         else  % Each frag. individually - HCD trigger - don't combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer group
%                 for j = 1:numfragmethods
%                     displaydataind_backup = displaydataind(i,j);
%                     displaydataind(i,j) = 0;
%                     if displaydataind_backup > 0
%                         if filterbyes
%                             tempenscore = allenscore(displaydataind_backup);
%                             thisfragmethod = allfragmode(displaydataind_backup);
%                             tempescutoff = filterbyescutoff(ismember(fragmodes,thisfragmethod));
%                             if tempenscore < tempescutoff
%                                 continue
%                             end
%                         end
%                         sgpseq = allsgpseq(displaydataind_backup);
%                         keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                             filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                         if ~keepthis
%                             continue
%                         end
%                         displaydataind(i,j) = displaydataind_backup;
%                     end
%                 end
%             end
%         end
%     else
%         if combineisomer  % Each frag. individually - Regular - combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer
%                 for j = 1:length(displaydataind{i})
%                     displaydataind_backup = displaydataind{i}(j);
%                     displaydataind{i}(j) = 0;
%                     if filterbyes
%                         tempenscore = allenscore(displaydataind_backup);
%                         thisfragmethod = allfragmode(displaydataind_backup);
%                         tempescutoff = filterbyescutoff(ismember(fragmodes,thisfragmethod));
%                         if tempenscore < tempescutoff
%                             continue
%                         end
%                     end
%                     sgpseq = allsgpseq(displaydataind_backup);
%                     keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                         filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                     if ~keepthis
%                         continue
%                     end
%                     displaydataind{i}(j) = displaydataind_backup;
%                 end
%             end
%         else  % Each frag. individually - Regular - don't combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer
%                 displaydataind_backup = displaydataind(i);
%                 displaydataind(i) = 0;
%                 if filterbyes
%                     tempenscore = allenscore(displaydataind_backup);
%                     thisfragmethod = allfragmode(displaydataind_backup);
%                     tempescutoff = filterbyescutoff(ismember(fragmodes,thisfragmethod));
%                     if tempenscore < tempescutoff
%                         continue
%                     end
%                 end
%                 sgpseq = allsgpseq(displaydataind_backup);
%                 keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                     filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                 if ~keepthis
%                     continue
%                 end
%                 displaydataind(i) = displaydataind_backup;
%             end
%         end
%     end
% else
%     % filter specified frag mode only
%     if scoreintdata.scoreoptions.isHCDtrigger
%         SASSO_method = scoreintdata.colnames;
%         methodind = find(ismember(SASSO_method,filterbyfragmode));
%         if isempty(methodind)
%             errordlg('Fragmentation method not found.','Fragmentation method not found.');
%         else
%             if combineisomer  % Specific frag mode - HCD trigger - combine isomer
%                 for i = 1:size(displaydataind,1)  % each isomer group
%                     for j = 1:length(displaydataind{i,1})  % each isomer
%                         displaydataind_backup = displaydataind{i,methodind}(j);
%                         displaydataind{i,methodind}(j) = 0;
%                         if displaydataind_backup > 0
%                             if filterbyes
%                                 tempenscore = allenscore(displaydataind_backup);
%                                 if tempenscore < filterbyescutoff
%                                     continue
%                                 end
%                             end
%                             sgpseq = allsgpseq(displaydataind_backup);
%                             keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                                 filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                             if ~keepthis
%                                 continue
%                             end
%                             displaydataind{i,methodind}(j) = displaydataind_backup;
%                         end
%                     end
%                 end
%             else  % Specific frag mode - HCD trigger - don't combine isomer
%                 for i = 1:size(displaydataind,1)
%                     displaydataind_backup = displaydataind(i,methodind);
%                     displaydataind(i,methodind) = 0;
%                     if displaydataind_backup > 0
%                         if filterbyes
%                             tempenscore = allenscore(displaydataind_backup);
%                             if tempenscore < filterbyescutoff
%                                 continue
%                             end
%                         end
%                         sgpseq = allsgpseq(displaydataind_backup);
%                         keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                             filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                         if ~keepthis
%                             continue
%                         end
%                         displaydataind(i,methodind) = displaydataind_backup;
%                     end
%                 end
%             end
%         end
%     else
%         if combineisomer  % Specific frag mode - Regular - combine isomer
%             for i = 1:size(displaydataind,1)  % each isomer group
%                 for j = 1:length(displaydataind{i})  % each isomer
%                     displaydataind_backup = displaydataind{i}(j);
%                     displaydataind{i}(j) = 0;
%                     if (displaydataind_backup > 0) && strcmpi(allfragmode{displaydataind_backup},filterbyfragmode)
%                         if filterbyes
%                             tempenscore = allenscore(displaydataind_backup);
%                             if tempenscore < filterbyescutoff
%                                 continue
%                             end
%                         end
%                         sgpseq = allsgpseq(displaydataind_backup);
%                         keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                             filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                         if ~keepthis
%                             continue
%                         end
%                         displaydataind{i}(j) = displaydataind_backup;
%                     end
%                 end
%             end
%         else  % Specific frag mode - Regular - don't combine isomer
%             for i = 1:size(displaydataind,1)
%                 displaydataind_backup = displaydataind(i);
%                 displaydataind(i) = 0;
%                 if (displaydataind_backup > 0) && strcmpi(allfragmode{displaydataind_backup},filterbyfragmode)
%                     if filterbyes
%                         tempenscore = allenscore(displaydataind_backup);
%                         if tempenscore < filterbyescutoff
%                             continue
%                         end
%                     end
%                     sgpseq = allsgpseq(displaydataind_backup);
%                     keepthis = checkquality(sgpseq,quality3,displaydataind_backup,...
%                         filterbyglypeponly,filterbyspecqual,filterbyprecion,filterbyuserdec);
%                     if ~keepthis
%                         continue
%                     end
%                     displaydataind(i) = displaydataind_backup;
%                 end
%             end
%         end
%     end
% end
% if scoreintdata.scoreoptions.isHCDtrigger  % Trigger experiment
%     numfragmethods = length(SASSO_method);
%     if combineisomer
%         for i = 1:size(displaydataind,1)
%             for j = 1:length(displaydataind{i,1})  % each isomer
%                 for k = 1:numfragmethods
%                     
%                 end
%             end
%         end
%     else
%         
%     end
% else
%     if combineisomer
%         
%     else
%         
%     end
% end
% [displaydata,~,~,~] = builddisptable(scoredata,displaydataind,displayoptions,usercustom);
% end