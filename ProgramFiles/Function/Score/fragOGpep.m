function [pepseq,pepcomp,pepmass,pepfrag_HCD,fragAAind_HCD,pepfrag_others,fragmode_others] = ...
    fragOGpep(protbatch,protindse,agpoptions)
thisprotbatch = protbatch(protindse(1):protindse(2));
pepfragcellnum = cellfun(@length,thisprotbatch);
pepseq = cell(sum(pepfragcellnum),1);
pepcomp = zeros(sum(pepfragcellnum),2);
pepmass = zeros(sum(pepfragcellnum),1);
fragmode_others = setdiff(agpoptions.fragmode,'HCD','stable');

counter = 1;
for ii = 1:length(thisprotbatch)
    tempprotbatch = thisprotbatch{ii};
    pepseq(counter:counter + length(tempprotbatch) - 1) = tempprotbatch(:);
    temppepcomp = ones(length(tempprotbatch),2) * ii + protindse(1) - 1;
    temppepcomp(:,2) = (1:length(tempprotbatch))';
    pepcomp(counter:counter + length(tempprotbatch) - 1,:) = temppepcomp;
    counter = counter + length(tempprotbatch);
end

pepseqind_ori = (1:length(pepseq))';
if agpoptions.doparallel
    parfor ii = 1:length(pepseq)
        pepmass(ii) = pepMW(pepseq{ii});
    end
    pepseqind_dist = distributed(pepseqind_ori);
    spmd
        [pepfrag_HCD_dist,fragAAind_HCD_dist,pepfrag_others_dist,pepseqindout_dist] = ...
            fragpepseq_og(pepseq,getLocalPart(pepseqind_dist),agpoptions);
    end
    pepseqind = [];
    pepfrag_HCD = {};
    fragAAind_HCD = {};
    pepfrag_others = {};
    for ii = 1:length(pepseqindout_dist)
        pepseqind = [pepseqind;pepseqindout_dist{ii}];
        pepfrag_HCD = [pepfrag_HCD;pepfrag_HCD_dist{ii}];
        fragAAind_HCD = [fragAAind_HCD;fragAAind_HCD_dist{ii}];
        pepfrag_others = [pepfrag_others;pepfrag_others_dist{ii}];
    end
    [~,ind] = sort(pepseqind);
    pepfrag_HCD = pepfrag_HCD(ind);
    fragAAind_HCD = fragAAind_HCD(ind);
    pepfrag_others = pepfrag_others(ind);
else
    pepmass = cellfun(@pepMW,pepseq);
    [pepfrag_HCD,fragAAind_HCD,pepfrag_others,~] = ...
        fragpepseq_og(pepseq,pepseqind_ori,agpoptions);
end
end

function [pepfrag_HCD,fragAAind_HCD,pepfrag_others,pepseqind] = ...
    fragpepseq_og(pepseq,pepseqind,agpoptions)
pepseqnow = pepseq(pepseqind);
fragmodes = agpoptions.fragmode;
pepfragnum = agpoptions.fragnum(:,1);
userpepiontyp = agpoptions.userpepiontyp;

HCDpepiontyp = userpepiontyp{strcmpi(fragmodes,'HCD')};
HCDpepnoniontyp = setdiff('bcyzi',HCDpepiontyp);
otherpepiontyp = userpepiontyp(~strcmpi(fragmodes,'HCD'));
otherfragmodes = setdiff(upper(fragmodes),'HCD','stable');
otherpepnoniontyp = cell(size(otherpepiontyp));
for ii = 1:length(otherpepiontyp)
    otherpepnoniontyp{ii} = setdiff('bcyzi',otherpepiontyp{ii});
end
[unipepfragnum,~,unipepfragnumind] = unique(pepfragnum);


pepfrag_HCD = cell(size(pepseqnow));
fragAAind_HCD = cell(size(pepseqnow));
pepfrag_others = cell(length(pepseqnow),length(otherpepiontyp));

for ii = 1:length(pepseqnow)
    if length(unipepfragnum) == 1
        temppepfrag = fragpep(pepseqnow{ii},unipepfragnum,'bcyzi','');
        temppepfrag_HCD = temppepfrag(~cellfun(@(x) any(ismember(x,HCDpepnoniontyp)),...
            {temppepfrag.type}));
        pepfrag_HCD{ii} = temppepfrag_HCD;
        tempfragAAind_HCD = cell(size(temppepfrag_HCD));
        for jj = 1:length(temppepfrag_HCD)
            tempfragAAind_HCD{jj} = temppepfrag_HCD(jj).unitindex{1};
        end
        fragAAind_HCD{ii} = tempfragAAind_HCD;
        for jj = 1:length(otherpepiontyp)
            pepfrag_others{ii,jj} = temppepfrag(~cellfun(@(x) any(ismember(x,otherpepnoniontyp{jj})),...
                {temppepfrag.type}));
        end
    elseif length(unipepfragnum) > 1
        tempunipepfrag = cell(size(unipepfragnum));
        for jj = 1:length(unipepfragnum)
            tempunipepfrag{jj} = fragpep(pepseqnow{ii},unipepfragnum(jj),'bcyzi','');
        end
        temppepfrag_HCD = tempunipepfrag{unipepfragnumind(strcmpi(fragmodes,'HCD'))};
        temppepfrag_HCD = temppepfrag_HCD(~cellfun(@(x) any(ismember(x,HCDpepnoniontyp)),...
            {temppepfrag_HCD.type}));
        pepfrag_HCD{ii} = temppepfrag_HCD;
        tempfragAAind_HCD = cell(size(temppepfrag_HCD));
        for jj = 1:length(temppepfrag_HCD)
            tempfragAAind_HCD{jj} = temppepfrag_HCD(jj).unitindex{1};
        end
        fragAAind_HCD{ii} = tempfragAAind_HCD;
        for jj = 1:length(otherpepiontyp)
            temppepfrag_other_jj = tempunipepfrag{unipepfragnumind(strcmpi(fragmodes,otherfragmodes{jj}))};
            pepfrag_others{ii,jj} = temppepfrag_other_jj(~cellfun(@(x) any(ismember(x,otherpepnoniontyp{jj})),...
                {temppepfrag_other_jj.type}));
        end
    end
end
end