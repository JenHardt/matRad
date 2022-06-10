function sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET,Pmod)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of Bragg peak degradation due to heterogeneities
%
% call
%   sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET)
%
% input
%   WET:        water equivalent thickness of heterogeneous structure, e.g. lung [mm]
%   Pmod:       modulation power [µm] (optional)
%
% output
%   sigmaSq:    sigma squared [mm^2] for degradation of Bragg peak
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% a = 1.43 [mm/sqrt(cm)] from Riccardos thesis, fit MC well,
% independent from proton energy, R80 stays at same position with degradation
% a = 1.60 [mm/sqrt(cm)] = 1.6/sqrt(10) [sqrt(mm)] fits the measurements better
%
% Pmod = a^2; 150-750 micrometer for swine lung from Witt et al.
% Pmod = 256 [µm] is equivalent to a = 1.6 [mm/sqrt(cm)]
% Pmod = 204.5 [µm] to a = 1.43 [mm/sqrt(cm)]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

if nargin < 2
    % "medium" modulation power
    % Pmod = 256; % [µm]

    % worst case modulation power
    Pmod = 800; % [µm]
end

% output used modulation power to console
matRad_cfg.dispInfo(['Modulation power set to Pmod = ' num2str(Pmod) ' µm.\n']);

% output sigma^2 in mm^2
sigmaSq = Pmod/1000 .* WET;
