function msdata = thinmsdata(msdata)
u = waitbar(0,'','Name','Thinning msdata...');
fields = fieldnames(msdata);
len = length(msdata.spectra);
for ll = 0:(len-1)
    if isempty(msdata.spectra{len-ll})
        hold_len = length(msdata.scannum);
        for ii = 1:length(fields)
            if length(msdata.(fields{ii})) == hold_len
                msdata.(fields{ii})(len-ll) = [];
            end
        end
    end
    waitbar((ll+1)/len,u,sprintf('%d/%d',(ll+1),len));
end
close(u)
end