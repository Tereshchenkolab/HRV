function [Pindex Gindex Eindex] = asymIndices(RR,N)

% Porta's and Guzik's asymmetry indices of a time series
% RR is the input time series
% N is the number of data points to be analysed

if nargin<2
  error('2 input arguments are required.')
end

if length(RR)<N
  error('Time points of input series are less than input data points, N.')
end

RR = RR(1:N);

% tau delayed time series
RRi = RR; RRi(end)=[];
RRitau = RR; RRitau(1)=[];

deltaRR = RRitau - RRi;

% Indices below and above the main diagonal of the plane
RRindices_up = find(deltaRR > 0); 
RRindices_down = find(deltaRR < 0);
RRindices_nonzero = find(deltaRR);

Pindex = [];
Pindex = (length(RRindices_down)/length(RRindices_nonzero))*100;

% Distance of time points from the main diagonal of the plane  
deltaRR_distance = (RRitau - RRi)/sqrt(2);

% dRRindices_up = find(deltaRR_distance > 0);
% dRRindices_down = find(deltaRR_distance < 0);

Gindex = [];
Gindex = (sum(deltaRR_distance(RRindices_up).^2)/sum(deltaRR_distance.^2))*100; 

Eindex = [];
Eindex = (sum(deltaRR.^3)/(sum(deltaRR.^2))^(3/2));
