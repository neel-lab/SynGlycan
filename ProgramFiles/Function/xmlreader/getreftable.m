function [subSASSO,followers_left,followerefs_left] = getreftable(starter,starteref,...
    followers,followerefs,columnseq,isreusable)
% starter: scan numbers, unique id's, this is the backbone of the table
% starterref: use these refs to find corresponding scans in "followers"
% follower: scan numbers, unique id's
% followerefs: refs of followers, search inside these to find which (Cell)
% follower to be associated with a member of starter (Cell too)
% columnseq: [starter, follower1, follower2, ...] put to which column

subSASSO = zeros(length(starter),length(followers)+1);

followerused = cell(size(followers));
for i = 1:length(followerused)
    followerused{i} = false(size(followers{i}));
end

for i = 1:length(starter)
    % triggered frag method with more scans first
    subSASSO(i,columnseq(1)) = starter(i);
    thisref = starteref(i);
    % the rest of the scans in this row
    for j = 1:length(followers)
        follower_scans = followers{j};
        follower_refscans = followerefs{j};
        availablefsn = follower_scans(~followerused{j});
        availablefrsn = follower_refscans(~followerused{j});
        identicalparentscan = availablefrsn == thisref;
        if any(identicalparentscan)
                tempind = find(identicalparentscan);
            subSASSO(i,columnseq(j+1)) = availablefsn(tempind(1));
            if ~isreusable(j)
                followerused{j}(follower_scans == availablefsn(tempind(1))) = true;
            end
        end
    end
end

followers_left = cell(size(followers));
followerefs_left = cell(size(followerefs));
for i = 1:length(followers)
    if any(~followerused{i})
        followers_left{i} = followers{i}(~followerused{i});
        followerefs_left{i} = followerefs{i}(~followerused{i});
    end
end
end

