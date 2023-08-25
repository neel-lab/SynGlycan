function decoymz = Iondecoy(input,type,orimz,oriAAind,origlyind,orimodind,mode)
% IONDECOY: calculate m/z for fragment decoys.
%
% Syntax:
% decoymz = Iondecoy(input,type,orimz,oriAAind,origlyind,orimodind,mode)
%
% Input:
% input: amino acid sequence or glycan part of BREAKGLYPEP result or the
% whole result.
% type: amino acid sequence permutation methods, available methods are:
%     "FIRST": swap first and last amino acids.
%     "FIRSTTWO": swap first and last, second and second to last amino
%         acids.
%     "FLIP": flip the whole sequence.
%     "RANDOM": randomize the whole sequence.
% orimz: m/z of all theoretical fragment ions.
% oriAAind: amino acid indices of fragments.
% origlyind: monosaccharide indices of fragments.
% orimodind: non-glycan PTM indices of fragments.
% mode: "glycan", "pep" or "glypep". 
%
% Output:
% decoymz: m/z value for each fragment of decoy
% glycan/peptide/glycopeptide.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function: 
% SEQPERM 
%

decoymz = orimz;
switch lower(mode)
    case 'pep'  % input is just AA sequence
        AAseq = input;
        pepmassdiff = zeros(1,length(AAseq));
        decoyseq = seqperm(AAseq,type);
        for ii = 1:length(pepmassdiff)
            pepmassdiff(ii) = AAMW(decoyseq(ii)) - AAMW(AAseq(ii));
        end
        for ii = 1:length(orimz)
            AAind = oriAAind{ii};
            decoymz(ii) = decoymz(ii) + sum(pepmassdiff(AAind));
        end
    case 'glycan'  % input is glyMat structure
        maxmwchange = 50;
        msmwchange = [];
        for ii = 1:length(input)
            msmwchange =  [msmwchange;...
                randfixedsum(input(ii).len,1,0,-maxmwchange,maxmwchange)];
        end
        for ii = 1:length(orimz)
            MSind = origlyind{ii}(:,1);
            decoymz(ii) = decoymz(ii) + sum(msmwchange(MSind));
        end
    case 'glypep'
        % input is 1 x 2 cell array, 1st is AA sequence, 2nd is glyMat
        pepdmz = Iondecoy(input{1}.pep,'RANDOM',orimz,oriAAind,origlyind,orimodind,'pep');
        glydmz = Iondecoy(input{2},'RANDOM',orimz,oriAAind,origlyind,orimodind,'glycan');
        decoymz = pepdmz + glydmz - decoymz;
end
end

function newseq = seqperm(seq,type)
% permutes AA seq number, e.g. [1,2,3,4] -> [2,4,1,3]
newseq = seq;
switch upper(type)
    case 'FIRST'
        newseq(1) = seq(end);
        newseq(end) = seq(1);
    case 'FIRSTTWO'
        newseq(1) = seq(end);
        newseq(2) = seq(end - 1);
        newseq(end) = seq(1);
        newseq(end-1) = seq(2);
    case 'FLIP'
        newseq = flip(seq);
    case 'RANDOM'
        newseq = newseq(randperm(length(seq)));
end
end