function varptm = getOSynglycans(varptmdata)
for ii = length(varptmdata)-1:-1:2
    varptm(ii).mod = varptmdata{ii};
    varptm(ii).aaresidue = 'S,T';
    varptm(ii).protpos = 0;
    varptm(ii).isglymod = true;
end
varptm(1).mod = '<o>';
varptm(1).aaresidue = 'M';
varptm(1).protpos = 0;
varptm(1).isglymod = false;
end

%{
function varptm = getOSynglycans(varptmfile)
varptm_mat = importdata(varptmfile);
OGlycans = varptm_mat.textdata;
for ii = length(OGlycans):-1:2
    varptm(ii).mod = OGlycans{ii};
    varptm(ii).aaresidue = 'S,T';
    varptm(ii).protpos = 0;
    varptm(ii).isglymod = true;
end
varptm(1).mod = '<o>';
varptm(1).aaresidue = 'M';
varptm(1).protpos = 0;
varptm(1).isglymod = false;
end
%}