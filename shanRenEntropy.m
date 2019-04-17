function [shannonEntr, renyiEntr] = shanRenEntropy(RRms,byteLevel,alpha)

% byteLevel = 4; %4ms window

% count the observations
observationHist = hist(RRms, byteLevel);

% convert to a probability
probXVal = observationHist ./ sum(observationHist);

% compute shannon entropy
shannonEntr = - sum( probXVal .* log2(probXVal) );

% compute renyi entropy
% alpha=4;
renyiEntr=(1/(1-alpha))*log2(sum(probXVal.^alpha));