function varptm = getNSynglycans(varptmdata)
for ii = length(varptmdata):-1:2
    varptm(ii).mod = varptmdata{ii};
    varptm(ii).aaresidue = 'N';
    varptm(ii).protpos = 0;
    varptm(ii).isglymod = true;
end
varptm(1).mod = '<o>';
varptm(1).aaresidue = 'M';
varptm(1).protpos = 0;
varptm(1).isglymod = false;
varptm(2)=[];
end

%{
function varptm = getNSynglycans(varptmfile)
varptm_mat = importdata(varptmfile);
NGlycans = varptm_mat.textdata;
for ii = size(NGlycans,1):-1:2
    varptm(ii).mod = NGlycans{ii,1};
    varptm(ii).aaresidue = NGlycans{ii,2};
    varptm(ii).protpos = 0;
    varptm(ii).isglymod = true;
end
varptm(1).mod = '<o>';
varptm(1).aaresidue = 'M';
varptm(1).protpos = 0;
varptm(1).isglymod = false;
end
%}