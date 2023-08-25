classdef Glycan
    % GLYCAN: An object storing the properties of single monosaccharides
    %    
    %   Syntax:
    %       monoMass = Glycan.glycanMSMap(mono1let) 
    %       monoFormula = Glycan.glycanformulaMap(mono1let) 
    %
    %   Input:
    %       mono1let: single-letter code for the monosaccharide
    %       
    %   Output;
    %       monoMass:  monosaccharide mass
    %       monoFormula: monosaccharide formula
    %    
    %   Examples:
    %       hexMass    = Glycan.glycanMSMap('h')
    %       hexFormula = Glycan.glycanformulaMap('h')
    %                      
    %   Note:                    Glycan table:
    %                          __________________
    %                            code | type
    %                          __________________
    %                              h  | Hexose
    %                              n  | HexNAc
    %                              s  | NeuAc
    %                              a  | aniline derivatized sialic acid
    %                              i  | ??
    %                              m  | iodoacetamide derivatized sialic acid 
    %                              g  | NeuGc
    %                              f  | Fuc
    %                              x  | Xyl
    %                              z  | SO3
    %                              p  | PO3
    %                              u  | KDN
    %                              k  | HexA
    %                              q  | GalNTGc
    %                              b  | 2AB derivatized sialic acid
    %                              d  | dimethylamine derivatized sialic acid
    %                              l  | sialolactone 
    %                              2e | initial n (N-glycans) without 2 molecules of water
    %                              e1 | initial n (N-glycans) without 1 molecule of water
    %                         _________________
    %
    %See also Anomeric
    
    %Author: Sriram Neelamegham
    %Date Lastly Updated: 1/28/2018
    
      properties (Constant)
             Hex         = struct('C',6,'H',10,'O',5);
             HexNAC      = struct('C',8,'H',13,'O',5,'N',1);
             Fuc         = struct('C',6,'H',10,'O',4);
             DeoHex      = struct('C',6,'H',10,'O',4);
             NeuAC       = struct('C',11,'H',17,'O',8,'N',1);
             AnaNeuAC    = struct('C',17,'H',22,'O',7,'N',2);  
             NeuGC       = struct('C',11,'H',17,'O',9,'N',1);        
             Pentose     = struct('C',5,'H',8,'O',4);                           
             Sulf        = struct('C',0,'H',0,'O',3,'S',1);  
             Phos        = struct('C',0,'H',1,'O',3,'P',1);               
             HexA        = struct('C',6,'H',8,'O',6);      
             KDN         = struct('C',9,'H',14,'O',8);
             Xyl         = struct('C',5,'H',8,'O',4);
             SO3         = struct('S',1,'O',3);
             PO3         = struct('P',1,'O',3,'H',1);
             HexNTGc     = struct('C',10,'H',16,'O',6,'N',2,'S',1);
             a23Sia      = struct('C',11,'H',15,'O',7,'N',1); % v
             a26Sia      = struct('C',12,'H',19,'O',8,'N',1); % w
             redEndForm        = struct('C',0,'H',2,'O',1,'N',0,'S',0,'P',0); 
             redEndMethylForm  = struct('C',2,'H',6,'O',1,'N',0,'S',0,'P',0);  %+2(ch2)
             HexMethyl         = struct('C',9,'H',16,'O',5,'N',0,'S',0,'P',0)  ;             %   +3(CH2)  3C6H
             HexNACMethyl      = struct('C',11,'H',19,'O',5,'N',1,'S',0,'P',0);     %   +3(CH2)  3C6H
             DeoHexMethyl      = struct('C',8,'H',14,'O',4,'N',0,'S',0,'P',0);    %  +2(CH2)  2C4H
             NeuACMethyl       = struct('C',16,'H',27,'O',8,'N',1,'S',0,'P',0);         %  +5(CH2)  5C10H
             
             
            glyfullname = {'Hex','HexNAc','NeuAc','AnaNeuAC','NeuGc','Fuc',...
                           'Xyl','SO3','PO3','KDN','HexA',...
                           'HexNTGc'};
                       
            gly1let     = {'h','n','s','a','i','m','g','f',...
                           'x','z','p','u','k',...
                           'q','b','d','v','w','l','e1','2e'};
                        
            glyformula = { struct('C',6,'H',10,'O',5),...
                           struct('C',8,'H',13,'O',5,'N',1),...
                           struct('C',11,'H',17,'O',8,'N',1),...
                           struct('C',17,'H',22,'O',7,'N',2),...  % <a>
                           struct('C',10,'H',15,'O',8,'N',1),...  % <i>
                           struct('C',13,'H',20,'O',9,'N',2),...  % <m>
                           struct('C',18,'H',23,'O',8,'N',3),...  % <b>
                           struct('C',11,'H',17,'O',9,'N',1),...
                           struct('C',6,'H',10,'O',4),...
                           struct('C',5,'H',8,'O',4),...
                           struct('S',1,'O',3),...
                           struct('P',1,'O',3,'H',1),...
                           struct('C',9,'H',14,'O',8),...
                           struct('C',6,'H',8,'O',6),...
                           struct('C',10,'H',16,'O',6,'N',2,'S',1),...
                           struct('C',13,'H',22,'O',7,'N',2),...    % <d>
                           struct('C',11,'H',15,'O',7,'N',1),...     % v
                           struct('C',12,'H',19,'O',8,'N',1),...     % w
                           struct('C',11,'H',15,'O',7,'N',1),...    % <l>
                           struct('C',8,'H',11,'O',4,'N',1),...    % <e1>
                           struct('C',8,'H',9,'O',3,'N',1)};    % <2e>
            
            glycanformulaMap = containers.Map({'h','n','s','a','i','m','g','f',...
                           'x','z','p','u','k',...
                            'q','b','d','v','w','l','e1','2e'},...
                          {struct('C',6,'H',10,'O',5),...
                           struct('C',8,'H',13,'O',5,'N',1),...
                           struct('C',11,'H',17,'O',8,'N',1),...
                           struct('C',17,'H',22,'O',7,'N',2),...
                           struct('C',10,'H',15,'O',8,'N',1),...  % <i>
                           struct('C',13,'H',20,'O',9,'N',2),...  % <m>
                           struct('C',18,'H',23,'O',8,'N',3),...  % <b>
                           struct('C',11,'H',17,'O',9,'N',1),...   
                           struct('C',6,'H',10,'O',4),...
                           struct('C',5,'H',8,'O',4),...
                           struct('S',1,'O',3),...
                           struct('P',1,'O',3,'H',1),...
                           struct('C',9,'H',14,'O',8),...
                           struct('C',6,'H',8,'O',6),...
                           struct('C',10,'H',16,'O',6,'N',2,'S',1),...
                           struct('C',13,'H',22,'O',7,'N',2),...  % <d>
                           struct('C',11,'H',15,'O',7,'N',1),...     % v
                           struct('C',12,'H',19,'O',8,'N',1),...     % w
                           struct('C',11,'H',15,'O',7,'N',1),...    % <l>
                           struct('C',8,'H',11,'O',4,'N',1),...    % <e1>
                           struct('C',8,'H',9,'O',3,'N',1)});    % <2e>
         
        glycanMSMap = containers.Map({'h','n','s','a','i','m','g','f',...
                           'x','z','p','u','k',...
                            'q','b','d','v','w','l','e1','2e','r','t'},...
           {162.0528235,203.0793724,291.0954164,366.14270,291.0954164,348.116883,307.0903311,146.0579089,...  % <i> is same as <s>
            132.0422588,79.95681459,79.96633093,176.0320881,250.0688675,...
            235.0514434,409.14851,318.14270,273.08484,305.11105,273.08485,185.06881,167.05824,306.10632,218.09027}); % GalNTGc=n+S

        glycanMSMeMap = containers.Map({'h','n','s','a','i','M','g','f',...
                           'x','z','p','u','k',...
                            'q','b','d','l','e1','2e'},...
           {204.0998,245.1263,361.1737,436.22098,347.1580508,432.210783,391.1842,174.0892,...   % note: <i>=<s>-14
            160.0736,79.95681459,79.96633093,218.0790,320.1472,...
            291.11402,507.25807,388.22095,315.13180,213.10011,181.07389});  % MW of GalNTGc is set to '0'. Remember to change it later on
        
         glycanMSAcetylMap = containers.Map({'h','n','s','a','i','m','g','f',...
                           'x','z','p','u','k',...
                            'q','b','d','e1'},...
           {288.0845,287.1005,417.1271,492.17440,417.1271,474.148578,475.1326, 230.0790,...   % <i> is same as <s>; should be edited later
            216.0634,79.95681459,79.96633093,260.0532,376.1006,...
            0,535.18021,486.18496,227.07937});  % MW of GlcNTGc is set to '0'. Remember to change it        
        
             %redEnd
             % redEnd =18.0105546;
             % redEndMethyl = 46.0419;
      end
end