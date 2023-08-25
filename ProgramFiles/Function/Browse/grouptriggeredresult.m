function tableindex = grouptriggeredresult(SASSO,protinds,scans,workingmode)
% GROUPTRIGGEREDRESULT: rearrange scoring results from triggered experiment
%
% Syntax:
% tableindex = grouptriggeredresult(SASSO,protinds,scans)
%
% Input:
% SASSO: m x n double. Each row represents a trigger - triggered
%     scan number combination, structure of the row is:
%     [T, T'd 1, T'd 2,..., PMS1, RT_PMS1]
%     Here:
%     T: "Trigger" scan number (HCD)
%     T'd 1, T'd 2,...: "Triggered" scan number (CID, ETHCD,...)
%     PMS1: scan number of parent MS1 scan
%     RT_PMS1: retention time of this parent MS1 scan.
% protinds: 1 x n cell array of strings. The compent serial numbers. Each
%     glycopeptide has a unique serial number. Format is [Protein, Peptide,
%     PTM1, PTM2,...].
% scans: 1 x n double. Scan number of all glycopeptides in the result.
%
% Output:
% tableindex: m x n double or m x n cell array of p x 1 double. Serial
%     number of identified glycopeptide to be placed in the final table. If
%     user choose not to group the isomers, this variable is double, meaning
%     the final table has m rows, n equals to the width of input "SASSO".
%     If user choose to combine isomers, this variable is cell array. p
%     equals to the number of isomers that are of this composition.
%
% Note:
% N/A
%
% Example:
% Set
%
% Children function:
% N/A
% See also:
% BUILDTRIGGERDISPTABLE
%

protinds = protinds(:);
scans = scans(:);
tblw = size(SASSO,2)-2;  % table width
switch upper(workingmode)
    case 'REGULAR'
        tableindex = [];
    case 'COMBINEISOMER'
        tableindex = {};
end

for ii = 1:size(SASSO,1)  % each group
    tempSASSOrow = SASSO(ii,1:tblw);
    tempresind = cell(1,tblw);
    protindref = cell(1,tblw);
    protindref2 = {};
    for jj = 1:tblw  % each frag method
        matchind = scans == tempSASSOrow(jj);
        tempresind{jj} = find(matchind);
        if any(matchind)
            protindref{jj} = [protindref{jj};protinds(matchind)];
            protindref2 = [protindref2;protinds(matchind)];
        end
    end
    if any(~cellfun(@isempty,protindref))  % found sgp
        uniprotindref = unique(protindref2);
        subtbl = zeros(length(uniprotindref),tblw);
        for jj = 1:tblw  % each fragmode
            tempprotindrefs = protindref{jj};
            for kk = 1:length(uniprotindref)  % each glypep
                scanispresentind = ismember(tempprotindrefs,uniprotindref{kk});
                if any(scanispresentind)
                    subtbl(kk,jj) = tempresind{jj}(scanispresentind);
                end
            end
        end
        switch upper(workingmode)
            case 'REGULAR'
                % INTENTIONALLY LEFT BLANK
            case 'COMBINEISOMER'
                subtbl_duplicate = subtbl;
                subtbl = cell(1,tblw);
                for jj = 1:tblw
                    subtbl{jj} = subtbl_duplicate(:,jj);
                end
        end
        tableindex = [tableindex;subtbl];
    end
end
keeprow = true(size(tableindex,1),1);
if isnumeric(tableindex)
    for ii = 1:size(tableindex,1)
        temptableindex_iszero = tableindex(ii,:) == 0;
        if any(temptableindex_iszero)
            keeprow(ii) = false;
        end
    end
elseif iscell(tableindex)
    for ii = 1:size(tableindex,1)
        temptableindex_iszero = false(1,length(tableindex{ii,1}));
        for jj = 1:tblw
            localtableindex = tableindex{ii,jj};
            localtableindex = reshape(localtableindex,1,[]);
            temptableindex_iszero = temptableindex_iszero | (localtableindex == 0);
        end
        if ~any(~temptableindex_iszero)
            keeprow(ii) = false;
        else
            for jj = 1:tblw
                tableindex{ii,jj} = tableindex{ii,jj}(~temptableindex_iszero);
            end
        end
    end
end
tableindex = tableindex(keeprow,:);
end