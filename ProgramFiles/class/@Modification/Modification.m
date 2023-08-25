classdef Modification
    %MODIFICATION: An object storing the properties of peptide modification
    %
    %  Syntax:
    %     modformula = Modification.formulaMap(mod1let);
    %
    %  Input:
    %     mod1let: single-letter for peptide modification type
    %
    %  Output:
    %     modformula: chemical formula for peptide modification
    %
    %  Example:
    %     modform = Modification.formulaMap('o');
    %
    %  Note:                    Modification table:
    %                   ___________________________________
    %                   single-letter | type
    %                   ___________________________________
    %                              o  | oxidation
    %                              s  | sulfation
    %                              i  | carbamidomethylation
    %                              c  | carboxymethylation
    %                              p  | Phosphorylation
    %                              a  | acetylation
    %                              b  | biotinylation
    %                              q  | Pyro-glu Q
    %                              e  | Pyro-glu E
    %                              f  | Formylation
    %                              n  | Sodiation
    %                              d  | Deamidation
    %                  ___________________________________
    %
    % See also: CHEMELE  GLYCAN  CHEMFORMULA  MODIFICATION  AMINOACID
    %     RESULTITEM
    
    % Author: Gang Liu
    % Date: 08/04/15
    % Last Updated: Kai Cheng
    % Date: 11/10/20
    
    properties (Constant)
        % Single letter code for modification
        mod1let = {'o','s','i',...
            'c','p','a',...
            'b','q','e',...
            'f','n','d'};
        
        % Full name of modification
        modfullname = {'oxidation','sulfation','carbamidomethylation',...
            'carboxymethylation','Phosphorylation','acetylation',...
            'biotinylation','Pyro-glu Q ',' Pyro-glu E',...
            'Formylation','Sodiation','Deamidation'};
        
        % Effect on amino acid mass. 1 for increase, 0 for loss
        modadd = {1,1,1,...
            1,1,1,...
            1,0,0,...
            1,1,1};
        
        % Formula of modifications. Sequence is identical to "mod1let" and
        %     "modfullname"
        modformula = {struct('O',1),...
            struct('S',1,'O',3),...
            struct('C',2,'H',3,'N',1,'O',1),...
            struct('C',2,'H',2,'O',2),...
            struct('P',1,'O',3,'H',1),...
            struct('C',2,'H',2,'O',1),...
            struct('C',10,'H',14,'N',2,'O',2,'S',1),...
            struct('N',1,'H',3),...
            struct('H',2,'O',1),...
            struct('C',1,'O',1),...
            struct('Na',1),...
            struct('O',1,'N',-1,'H',-1)};
        
        % Retrieve modification composition using single letter codes
        formulaMap = containers.Map({'o','s','i',...
            'c','p','a',...
            'b','q','e',...
            'f','n','d'},...
            {struct('O',1),...
            struct('S',1,'O',3),...
            struct('C',2,'H',3,'N',1,'O',1),...
            struct('C',2,'H',2,'O',2),...
            struct('P',1,'O',3,'H',1),...
            struct('C',2,'H',2,'O',1),...
            struct('C',10,'H',14,'N',2,'O',2,'S',1),...
            struct('N',1,'H',3),...
            struct('H',2,'O',1),...
            struct('C',1,'O',1),...
            struct('Na',1),...
            struct('O',1,'N',-1,'H',-1)});
        
        % Retrieve effect on amino acid mass using single letter codes
        mwEffectMap = containers.Map({'o','s','i',...
            'c','p','a',...
            'b','q','e',...
            'f','n','d'},...
            {1,1,1,...
            1,1,1,...
            1,0,0,...
            1,1,1});
        
        % RESERVED - Isobaric tags
        isotagname = {'iTRAQ4_114';'iTRAQ4_115';'iTRAQ4_116';'iTRAQ4_117';...
            'iTRAQ8_113';'iTRAQ8_114';'iTRAQ8_115';'iTRAQ8_116';...
            'iTRAQ8_117';'iTRAQ8_118';'iTRAQ8_119';'iTRAQ8_121';...   % 8-plex iTRAQ
            'TMT0_126';...  % 0-plex TMT
            'TMT2_126';'TMT2_127';...  % 2-plex TMT
            'TMT6_126';'TMT6_127';'TMT6_128';'TMT6_129';...
            'TMT6_130';'TMT6_131';'TMT10_126';'TMT10_127N';...
            'TMT10_127C';'TMT10_128N';'TMT10_128C';'TMT10_129N';...
            'TMT10_129C';'TMT10_130N';'TMT10_130C';'TMT10_131';...
            'TMT10_131C'};  % 6/10-plex TMT
    end
    
    methods(Static)
        % Get single letter code for modification
        function modchar = getmodchar
            modchar = char(Modification.mod1let)';
        end
        
        % Create expression for regular expression matching
        function aaexpr  = getmodcharexpr
            aaexpr =['[^' Modification.getmodchar];
            aaexpr =[aaexpr ']'];
        end
        
        % Check if input modification single letter code is available for use
        function ismodchar = ismodtring(stringinput)
            findnonmodchar = regexp(stringinput,Modification.getmodcharexpr, 'once');
            ismodchar      = isempty(findnonmodchar);
        end
    end
    
end
