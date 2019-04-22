
function xf = filtbyfft(x, fs, passband)
%FILTBYFFT: filter by FFT
% x: signal vector, fs: sampling frequency,  passband: passband interval
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

if nargin<3
  error('3 input arguments are required. ')
end

%if fs<=2*fcut
%  error('Cut-off frequency must be smaller than half of sampling frequency.')
%end

if length(passband)~=2
  error('Passband must be a vector of two components specifying the passband interval.')
end

passband(1) = max(passband(1), 0);
passband(2) = min(passband(2), 0.5*fs);

N = length(x);
y = fft(x);

lowicut = round(passband(1)*N/fs);
lowmirror = N-lowicut+2;

highicut =  round(passband(2)*N/fs);
highmirror = N-highicut+2;

y([1:(lowicut-1) (lowmirror+1):end])=0;
y((highicut+1):(highmirror-1))=0;


xf = ifft(y);