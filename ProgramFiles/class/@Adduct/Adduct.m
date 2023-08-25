classdef Adduct
    % Adduct: An object storing the properties of the adduct entity
    %    
    %   Syntax:
    %       adductMass = Adduct.mass(adductCell) 
    %       adductFormula = Adduct.formulaMap(adductCell) 
    %
    %   Input:
    %       adductCell: Cell describing the adduct
    %       
    %   Output;
    %       adductMass:  Adduct mass
    %       adductFormula: Adduct formula
    %    
    %   Examples:
    %       adduct{1}.add='Na';
    %       adduct{1}.count=2;    
    %       Mass    = Adduct.mass(adduct(1}.add)
    %       Formula = Adduct.formulaMap(adduct{1}.add)
    %                      
    %   Note:                 Anomeric end table:
    %                          __________________
    %                            code | type
    %                          __________________
    %                              H  | hydrogen
    %                              Li | Lithium
    %                              Na | Sodium
    %                              K  | Pottasium
    %                              Cl | Chlorine
    %                              Br | Bromine
    %                             NH4 | Ammonium
    %                          _________________
    %
    %See also Anomeric, Glycan.
    
    %Author: Sriram Neelamegham
    %Date Lastly Updated: 4/23/2021
    
      properties (Constant)
             H         = struct('H',1);
             Li        = struct('Li',1);
             Na         = struct('Na',1);
             K        = struct('K',1);
             Cl        = struct('Cl',1);
             Br        = struct('Br',1);
             NH4       = struct('C',0,'N',1,'H',4);
             
             adductfullname = {'H','Li','Na','K','Cl','Br','NH4'};
             adduct1let     = {'H','L','N','K','C','B','A'};
             adductformula =  {struct('H',1)...
                           struct('Li',1)...
                           struct('Na',1)...
                           struct('K',1)...
                           struct('Cl',1)...
                           struct('Br',1)...
                           struct('N',1,'H',4)};  
              adductcharge = [1,1,1,1,-1,-1,1];
            
              formulaMap = containers.Map({'H','Li','Na','K','Cl','Br','NH4'},...
                           {struct('H',1)...
                           struct('Li',1)...
                           struct('Na',1)...
                           struct('K',1)...
                           struct('Cl',1)...
                           struct('Br',1)...
                           struct('N',1,'H',4)});
           
               mass = containers.Map({'H','Li','Na','K','Cl','Br','NH4'},...
                {1.00727646677,7.0160041,22.98976966,38.9637069,34.96885271,78.9183379,18.038466}); 
            
               charge = containers.Map({'H','Li','Na','K','Cl','Br','NH4'},...
                {1,1,1,1,-1,-1,1}); 
      end
end