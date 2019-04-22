function [picind, s] = rpeak(s0, fs)
% RPEAKS: detection of ECG R-peaks by parabolic fitting
% [R, SF] = RPEAKS(S,FS), S is the signal vector, FS is the sampling frequency in HZ.
% R is the indices of the R-peaks, SF is the filtered signal
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

s = filtbyfft(s0-mean(s0),fs, [0.5 300]);
N=length(s);

fratio = fs/300; % Frequence ratio w.r.t. 250hz.

hwin=ceil(5*fratio); % half window width for parabolic fitting
weights=1:(-1/hwin):(1/hwin); % parabolic fitting weights in half window
dsw = sqrt([fliplr(weights), 1, weights])'; % square root of double side weights

ps2 = zeros(size(s));
ps0 = zeros(size(s));

xw = [-[(-hwin:hwin).^2]', ones(2*hwin+1,1)];
xw = dsw(:,[1 1]) .* xw;

invxw = inv(xw'*xw)*xw';
for k=(hwin+1):(N-hwin)
  yw = s((k-hwin):(k+hwin));
  yw = dsw .* yw;
  
  %th = xw\yw;
  th = invxw*yw;
  
  ps2(k) = th(1);
  ps0(k) = th(2);
end
pth = ps2.*ps0;
% END of indicator signal computation

% R location
%---------------
PTc = zeros(N,1);
kpt = 0;

pth0 = pth; % Store original pth for missed peaks recovery

%sortedpth = sort(pth);
% Added floor to convert to integer because it's used as an index and was
% throwing errors otherwise
srtlen = floor(min(10000*fratio, length(pth)));
sortedpth = sort(pth(1:srtlen));

%seuil = sortedpth(round(N*0.975));
%%seuil = sortedpth(round(srtlen*0.975));

%Due to small amplitudes in R peaks (i.e. 13 MI data, Lead I, at 13 sec)
seuil = sortedpth(round(srtlen*0.975));  

pth(find(pth<seuil)) = 0;


cnth = round(0.36*fs);
debut = 1;
fin = 0;
for k=(1+cnth):N-cnth
  if pth(k) & ~any(pth((k-cnth):(k-1)))
    debut = k;
  end
  if pth(k) & ~any(pth((k+1):(k+cnth)))
    fin = k;
    [dum, ind] = max(pth(debut:fin));
    ind = ind + debut - 1;
    kpt = kpt+1;
    PTc(kpt) = ind;
  end
end
picind = PTc(1:kpt);

if length(picind)>1
    % Look for missed peaks
    srtlen = min(10000*fratio, kpt);
    dpic = diff(picind(1:srtlen));
    sorteddpic = sort(dpic);
    typicnum = median(sorteddpic(1:ceil(srtlen/3)));

    %Treat possible beginning peaks
    beginexind = BeginPeakRecovery(picind(1), pth0, typicnum, seuil);
    exind = EndPeakRecovery(picind(end), pth0, typicnum, seuil,N);
    %exind=[];
    picind = [sort(beginexind); picind; sort(exind)];

    % These 3 lines have to be redone since now picind may have been changed.
    dpic = diff(picind(1:srtlen));
    sorteddpic = sort(dpic);
    typicnum = median(sorteddpic(1:ceil(srtlen/3)));

    indmiss = find(dpic>typicnum*1.5);
    if ~isempty(indmiss)
      misscase = length(indmiss);
      extraind = cell(misscase,1);
      for km=1:length(indmiss)
        extraind{km} = PeakRecovery(picind(indmiss(km)), picind(indmiss(km)+1), pth0, typicnum);
      end
      vecextraind = cell2mat(extraind);
      keptind = find(pth0(vecextraind)>0.3*seuil);
      picind = sort([picind; vecextraind(keptind)]);
    end
end

%kpt = length(picind);


%=================================================================
function exind = PeakRecovery(ind1, ind2, pth, typicnum)
% Recover missed peaks between two detected peaks

hnum = ceil(0.5*typicnum);
[dum, maxind] = max(pth((ind1+hnum):(ind2-hnum)));
maxindglobal = maxind + ind1+hnum - 1;
exind =  maxindglobal;

if (maxindglobal-ind1)>typicnum*1.5
  exind = [PeakRecovery(ind1, maxindglobal, pth, typicnum); exind];
end

if (ind2-maxindglobal)>typicnum*1.5
  exind = [exind; PeakRecovery(maxindglobal, ind2, pth, typicnum)];
end 

%-----------------------------------
function exind = BeginPeakRecovery(ind1, pth0, typicnum, seuil)
% Recover missed peaks at the beginning of the signal

if ind1>typicnum*0.7 
  [maxv, maxi] = max(pth0(1:ceil(ind1-typicnum*0.3)));
  if maxv>seuil*0.8
    exind = [maxi; BeginPeakRecovery(maxi, pth0, typicnum, seuil)];
  else
    exind = [];
  end
else
  exind = [];
end
%%

function exind = EndPeakRecovery(ind2, pth0, typicnum, seuil,N)
% Recover missed peaks at the beginning of the signal

if N-ind2>typicnum %0.7
  [maxv, maxi] = max(pth0(ceil(ind2+typicnum*0.3):end));
% if maxv>seuil*0.6
%     exind = [maxi; EndPeakRecovery(maxi, pth0, typicnum, seuil,N)];
%   else
%     exind = [];
%   end
exind=ceil(ind2+typicnum*0.3)+maxi-1;
% else
%   exind = [];
% end
else
    exind =[];
end


%==========================================================================
