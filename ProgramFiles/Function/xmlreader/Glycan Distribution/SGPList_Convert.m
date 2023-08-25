function SGPList = SGPList_Convert(SGPList)
mass_h = 180.156;
mass_n = 221.209;
mass_s = 268.218;
mass_f = 164.157;
for ii = 1:length(SGPList)
    SGPList{ii,2} = count(SGPList{ii,1},'h');
    SGPList{ii,3} = count(SGPList{ii,1},'n');
    SGPList{ii,4} = count(SGPList{ii,1},'s');
    SGPList{ii,5} = count(SGPList{ii,1},'f');
    SGPList{ii,6} = (SGPList{ii,2}+SGPList{ii,3}+SGPList{ii,4}+SGPList{ii,5});
    SGPList{ii,7} = (SGPList{ii,2}*mass_h)+(SGPList{ii,3}*mass_n)+(SGPList{ii,4}*mass_s)+(SGPList{ii,5}*mass_f);
end
end