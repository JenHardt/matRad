function samples = matRad_sampleBino(n,p,numOfSamples,continuous)
% matRad function for binomial sampling
%
% call
%   X = matRad_sampleBino(n,p,numOfSamples)
%
% input
%   n:              number of independent experiments
%   p:              probability (between 0 and 1)
%   numOfSamples:   number of samples for output
%
% output
%   X:              binomial samples
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

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% handle inputs
if any(~mod(n,1) == 0) && ~continuous
    matRad_cfg.dispError('n has to be integer for discrete distribution');
elseif ~(all(p <= 1) && all(p >= 0))
    matRad_cfg.dispError('p must be between 0 and 1');
elseif ~isscalar(p) && ~(numOfSamples == numel(p))
    matRad_cfg.dispError('p array must have numOfSamples entries')
elseif ~isscalar(n) && ~(numOfSamples == numel(n))
    matRad_cfg.dispError('n array must have numOfSamples entries')
end

% save time when testing on homogeneous phantom or when performing density
% override
if isscalar(unique(n)) && isscalar(unique(p))
    n = unique(n);
    p = unique(p);
end

if continuous
    % calculate beta distribution parameters alpha & beta using method of moments with sample mean x=p and 
    % https://en.wikipedia.org/wiki/Beta_distribution#Method_of_moments
    a = p .* (n-1);
    b = (1-p) .* (n-1);

%     a = p.*n;
%     b = (1-p).*n;

    % sample from continuous beta distribution using "numOfSamples" random numbers
    samples = betaincinv(rand([numOfSamples,1]),a,b);
else
    % sample discrete binomial distribution
    if isscalar(n)
        samples = sum(rand([numOfSamples,n]) < p, 2);
    else
        samples = zeros(numOfSamples,1);
        for i = 1:numel(n)
            samples(i) = sum(rand([1,n(i)]) < p(i), 2);
        end
    end
    % need normalization here to only get values over [0,1]
    samples = samples ./ n;
end

end
