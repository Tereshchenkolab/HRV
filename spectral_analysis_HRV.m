function [LF,HF,VLF,TP,LFHFratio] = spectral_analysis_HRV(RR,Fs,type)
%spectral_analysis_HRV Spectral analysis of HRV.
%   The function uses FFT to compute the spectral density function of the
%   interpolated RR tachogram.
%
% Input:
%   RR: RR intervals in seconds.
%   Fs: Sampling frequency.
%   type: interpolation type. Look up interp1 function of Matlab for
%   accepted types (default: 'spline').
%
% Output:
%   LF: Low frequency Power
%   HF: High frequency Power
%   VLF: Very low frequency Power
%   TP: Total Power
%   LFHFratio: LF/HF
     
    if nargin<2 || isempty(Fs)
        error('wrong number or types of arguments');
    end 
    
    if nargin<3
        type = 'spline';
    end
    
    RR = RR(:);
    
    switch type
        case 'none'
            RR_rsmp = RR;
        otherwise
            if sum(isnan(RR))==0 && length(RR)>1
                ANN = cumsum(RR)-RR(1);
                % use interp1 methods for resampling
                RR_rsmp = interp1(ANN,RR,0:1/Fs:ANN(end),type);
            else
                RR_rsmp = [];
            end
    end
    
    % FFT
    L = length(RR_rsmp); 
    
    if L==0 
        LFHFratio = NaN;
        VLF = NaN;
        LF = NaN;
        HF = NaN;
        TP = NaN;
    else
        NFFT = 2^nextpow2(L);
        Y = fft(zscore(RR_rsmp),NFFT)/L;
        %Y = fft(nanzscore(RR_rsmp),NFFT)/L;
        f = Fs/2*linspace(0,1,NFFT/2+1);  

        YY = 2*abs(Y(1:NFFT/2+1));

        VLF = sum(YY(f<=0.04));
        LF = sum(YY(f<=0.15))-VLF;  
        HF = sum(YY(f<=0.4))-VLF-LF;
        TP = sum(YY(f<=0.4));
  
        LFHFratio = LF/HF; 
    end
end
