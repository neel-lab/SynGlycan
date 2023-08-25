classdef Aminoacid
    % AMINOACID: An Object storing the properties of amino acids
    %
    %   Syntax:
    %        aa3letlist = Aminoacid.aa3let;
    %        aa1letlist = Aminoacid.aa1let;
    %        aaformula = Aminoacid.formulaMap(aa1let);
    %
    %   Input:
    %        aa1let: single-letter code representing single amino acid
    %            (Upper case only)
    %
    %   Output:
    %        aa1letlist: single-letter code for all amino acids
    %        aa3letlist: 3-letter code for all amino acids
    %        aaformula: chemical formula for single amino acid
    %
    %   Example:
    %       Example 1:
    %          aa3letlist = Aminoacid.aa3let
    %
    %       Example 2:
    %          aa1letlist = Aminoacid.aa1let
    %
    %       Example 3:
    %          alaformula = Aminoacid.formulaMap('A');
    %
    %   Note:                    Amino Acid table:
    %                          __________________
    %                            code | type
    %                          __________________
    %                              A  |  Ala
    %                              R  |  Arg
    %                              N  |  Asn
    %                              D  |  Asp
    %                              C  |  Cys
    %                              Q  |  Gln
    %                              E  |  Glu
    %                              G  |  Gly
    %                              H  |  His
    %                              I  |  Ile
    %                              L  |  Leu
    %                              K  |  Lys
    %                              M  |  Met
    %                              F  |  Phe
    %                              P  |  Pro
    %                              S  |  Ser
    %                              T  |  Thr
    %                              W  |  Trp
    %                              Y  |  Tyr
    %                              V  |  Val
    %                              B  |  Asx
    %                              Z  |  Glx
    %                              X  |  Unk
    %                              U  |  Sec
    %                          _________________
    %
    %See also Chemele, Glycan, Chemformula, Protease, Modification, mzXML
    %
    % Author: Gang Liu
    % Date: 1/19/2015
    
    properties (Constant)
        % AA3LET: 3 letter code for amino acids
        aa3let = {'Ala','Arg',...  % 1
            'Asn','Asp',...  % 2
            'Cys','Gln',...  % 3
            'Glu','Gly',...  % 4
            'His','Ile',...  % 5
            'Leu','Lys',...  % 6
            'Met','Phe',...  % 7
            'Pro','Ser',...  % 8
            'Thr','Trp',...  % 9
            'Tyr','Val',...  % 10
            'Unk','Sec'};...  % 12
        
        % AA1LET: single letter code for amino acids
        aa1let = {'A','R',...  % 1
            'N','D',...  % 2
            'C','Q',...  % 3
            'E','G',...  % 4
            'H','I',...  % 5
            'L','K',...  % 6
            'M','F',...  % 7
            'P','S',...  % 8
            'T','W',...  % 9
            'Y','V',...  % 10
            'X','U'};...  % 12
            
        % AAFORMULA: composition of amino acids, sequence is identical to
        %     "AA3LET" and "AA1LET"
        aaformula = {struct('C',3,'H',7,'N',1,'O',2),struct('C',6,'H',14,'N',4,'O',2),...  % 1
            struct('C',4,'H',8,'N',2,'O',3),struct('C',4,'H',7,'N',1,'O',4),...  % 2
            struct('C',3,'H',7,'N',1,'O',2,'S',1),struct('C',5,'H',10,'N',2,'O',3),...  % 3
            struct('C',5,'H',9,'N',1,'O',4),struct('C',2,'H',5,'N',1,'O',2),...  % 4
            struct('C',6,'H',9,'N',3,'O',2),struct('C',6,'H',13,'N',1,'O',2),...  % 5
            struct('C',6,'H',13,'N',1,'O',2),struct('C',6,'H',14,'N',2,'O',2),...  % 6
            struct('C',5,'H',11,'N',1,'O',2,'S',1),struct('C',9,'H',11,'N',1,'O',2),...  % 7
            struct('C',5,'H',9,'N',1,'O',2),struct('C',3,'H',7,'N',1,'O',3),...  % 8
            struct('C',4,'H',9,'N',1,'O',3),struct('C',11,'H',12,'N',2,'O',2),...  % 9
            struct('C',9,'H',11,'N',1,'O',3),struct('C',5,'H',11,'N',1,'O',2),...  % 10
            struct('C',0,'H',0,'N',0,'O',0,'S',0),struct('C',3,'H',7,'N',1,'O',2,'Se',1)};  % 12
        
        % FORMULAMAP: retrieve amino acid composition with single letter
        %     code input
        formulaMap = containers.Map({'A','R',...  % 1
            'N','D',...  % 2
            'C','Q',...  % 3
            'E','G',...  % 4
            'H','I',...  % 5
            'L','K',...  % 6
            'M','F',...  % 7
            'P','S',...  % 8
            'T','W',...  % 9
            'Y','V',...  % 10
            'X','U'},...  % 12
            {struct('C',3,'H',7,'N',1,'O',2),struct('C',6,'H',14,'N',4,'O',2),...  % 1
            struct('C',4,'H',8,'N',2,'O',3),struct('C',4,'H',7,'N',1,'O',4),...  % 2
            struct('C',3,'H',7,'N',1,'O',2,'S',1),struct('C',5,'H',10,'N',2,'O',3),...  % 3
            struct('C',5,'H',9,'N',1,'O',4),struct('C',2,'H',5,'N',1,'O',2),...  % 4
            struct('C',6,'H',9,'N',3,'O',2),struct('C',6,'H',13,'N',1,'O',2),...  % 5
            struct('C',6,'H',13,'N',1,'O',2),struct('C',6,'H',14,'N',2,'O',2),...  % 6
            struct('C',5,'H',11,'N',1,'O',2,'S',1),struct('C',9,'H',11,'N',1,'O',2),...  % 7
            struct('C',5,'H',9,'N',1,'O',2),struct('C',3,'H',7,'N',1,'O',3),...  % 8
            struct('C',4,'H',9,'N',1,'O',3),struct('C',11,'H',12,'N',2,'O',2),...  % 9
            struct('C',9,'H',11,'N',1,'O',3),struct('C',5,'H',11,'N',1,'O',2),...  % 10
            struct('C',0,'H',0,'N',0,'O',0,'S',0),struct('C',3,'H',7,'N',1,'O',2,'Se',1)});...  % 12
    end

methods(Static)    
    % Get single letter amino acid codes
    function aachar = getaachar
        aachar = char(Aminoacid.aa1let)';
    end
    
    % Create expression for regular expression matching
    function aaexpr = getaacharexpr
        aaexpr =['[^' Aminoacid.getaachar];
        aaexpr =[aaexpr ']'];
    end
    
    % List of available amino acid single letter codes
    function aaexpr = getaa1letcharexpr
        aaexpr ='[ARNDCQEGHILKMFPSTWYVUO]';
    end
    
    % Check if input amino acid single letter code is available for use
    function isaachar = isaastring(stringinput)
        findnonaachar = regexp(stringinput,Aminoacid.getaacharexpr, 'once');
        isaachar      = isempty(findnonaachar);
    end
end
end

