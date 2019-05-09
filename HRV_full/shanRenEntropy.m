function [shannonEntr, renyiEntr] = shanRenEntropy(RR,byteLevel,alpha)
% shanRenEntropy: This function calculates Shannon Entropy and Renyi Entropy.
%
%   RR       : RR intervals in milliseconds.
%   byteLevel: bins (class width) for calculating histogram of RR
%   alpha    : weighting coefficient 
%

if nargin<3 || isempty(alpha)
    alpha=4; %weighs larger probabilities more than lower coefficients
% alpha = 3 or 4 for Renyi entropy improved accuracy 
% in separating congestive heart failure from normal sinus rhythm.
% Conforth et al., CinC 2016
end 

if nargin<2 || isempty(byteLevel)
    byteLevel = 4; %4ms window / width
end 

% count the observations
observationHist = hist(RR, byteLevel);

% convert to a probability
probXVal = observationHist ./ sum(observationHist);

% compute shannon entropy
shannonEntr = - sum( probXVal .* log2(probXVal) );

% compute renyi entropy from weight probability distribution within tachogram
renyiEntr=(1/(1-alpha))*log2(sum(probXVal.^alpha));
