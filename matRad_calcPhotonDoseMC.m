function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,calcDoseDirect)
% matRad Monte Carlo photon dose calculation wrapper
%   Will call the appropriate subfunction for the respective
%   MC dose-calculation engine
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst)
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nHistories:                 number of histories per beamlet
%   calcDoseDirect:             (optional) binary switch to enable forward 
%                               dose calcualtion (default: false)
% output
%   dij:                        matRad dij struct
%
% References
%
% -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

if nargin < 6
    calcDoseDirect = false;
end

% load appropriate config from pln or from class
pln = matRad_cfg.getDefaultClass(pln,'propMC');

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,'propDoseCalc');

switch pln.propMC.engine
    case 'ompMC'
        dij = matRad_calcPhotonDoseOmpMC(ct,stf,pln,cst,nHistories,calcDoseDirect);
    case 'TOPAS'
        dij = matRad_calcPhotonDoseMCtopas(ct,stf,pln,cst,calcDoseDirect);
    otherwise
        matRad_cfg.dispError('MC engine %s not known/supported!',engine);
end

end
