function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin2(ecg,fs)

%% function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs)
% Complete implementation of Pan-Tompkins algorithm

%% Inputs
% ecg : raw ecg vector signal 1d signal
% fs : sampling frequency e.g. 200Hz, 400Hz and etc

%% Outputs
% qrs_amp_raw : amplitude of R waves amplitudes
% qrs_i_raw : index of R waves

%% Method :

%% PreProcessing
% 1) Signal is preprocessed with the same filtering setups introduced in Pan
% tompkins paper (a combination of low pass and high pass filter 5-15 Hz)
% to get rid of the baseline wander and muscle noise. 

% 2) The filtered signal
% is derivated using a derivating filter to high light the QRS complex.

% 3) Signal is squared.4)Signal is averaged with a moving window to get rid
% of noise (0.150 seconds length).

% 5) depending on the sampling frequency of your signal the filtering
% options are changed to best match the characteristics of your ecg signal

% 6) Unlike the other implementations in this implementation the desicion
% rule of the Pan tompkins is implemented completely.

%% Decision Rule 
% At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
% waveform at the output of the MWI . The determination as to whether this pulse
% corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
% performed with an adaptive thresholding operation and other decision
% rules outlined below;

% a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
% samples at the location of the MWI maxima. This is done in order to localize the QRS
% complex to a single instant of time. The w[k] weighting is the maxima value.

% b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
% two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
% 2 second training phase) that continuously adapt to changing ECG signal quality. The
% first pass through y[n] uses these thresholds to classify the each non-zero sample
% (CURRENTPEAK) as either signal or noise:
% If CURRENTPEAK > THR_SIG, that location is identified as a �QRS complex
% candidate� and the signal level (SIG_LEV) is updated:
% SIG _ LEV = 0.125 �CURRENTPEAK + 0.875� SIG _ LEV

% If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
% �noise peak� and the noise level (NOISE_LEV) is updated:
% NOISE _ LEV = 0.125�CURRENTPEAK + 0.875� NOISE _ LEV
% Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
% respectively) at that point in the ECG, the thresholds are adjusted as follows:
% THR _ SIG = NOISE _ LEV + 0.25 � (SIG _ LEV ? NOISE _ LEV )
% THR _ NOISE = 0.5� (THR _ SIG)
% These adjustments lower the threshold gradually in signal segments that are deemed to
% be of poorer quality.


% c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
% CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
% complex. If however, an unreasonably long period has expired without an abovethreshold
% peak, the algorithm will assume a QRS has been missed and perform a
% searchback. This limits the number of false negatives. The minimum time used to trigger
% a searchback is 1.66 times the current R peak to R peak time period (called the RR
% interval). This value has a physiological origin - the time value between adjacent
% heartbeats cannot change more quickly than this. The missed QRS complex is assumed
% to occur at the location of the highest peak in the interval that lies between THR_SIG and
% THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is 
% calculated as an average of the last eight QRS locations in order to adapt to changing heart 
% rate and the second RR interval mean is the mean 
% of the most regular RR intervals . The threshold is lowered if the heart rate is not regular 
% to improve detection.

% d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
% impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
% detected one. This constraint is a physiological one � due to the refractory period during
% which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
% candidates are generated, the algorithm eliminates such physically impossible events,
% thereby reducing false positives.

% e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
% refractory period but within 360ms of the previous QRS, the algorithm determines
% whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
% T wave. This decision is based on the mean slope of the waveform at that position. A slope of
% less than one half that of the previous QRS complex is consistent with the slower
% changing behaviour of a T wave � otherwise, it becomes a QRS detection.
% Extra concept : beside the points mentioned in the paper, this code also
% checks if the occured peak which is less than 360 msec latency has also a
% latency less than 0,5*mean_RR if yes this is counted as noise

% f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
% checked with the help of the output of the bandpass signal to improve
% detection and find the original index of the real R waves on the raw ecg
% signal

%% References :

%[1]PAN.J, TOMPKINS. W.J,"A Real-Time QRS Detection Algorithm" IEEE
%TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985.

%% Author : Hooman Sedghamiz
% Linkoping university 
% email : hoose792@student.liu.se
% hooman.sedghamiz@medel.com

% Any direct or indirect use of this code should be referenced 
% Copyright march 2014
%%
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end


ecg = ecg(:); % vectorize

%% Initialize
qrs_c =[]; %amplitude of R
qrs_i =[]; %index
SIG_LEV = 0; 
nois_c =[];
nois_i =[];
NOISE_LEV = 0;
delay = 0;
skip = 0; % becomes one when a T wave is detected
not_nois = 0; % it is not noise when not_nois = 1
selected_RR =[]; % Selected RR intervals
m_selected_RR = 0;
mean_RR = 0;
qrs_i_raw =[];
qrs_amp_raw=[];
SIG_LEV= 0;
ser_back = 0; 
SIG_LEV1 = 0; % Signal level in Bandpassed filter
NOISE_LEV1 = 0; % Noise level in Bandpassed filter
test_m = 0;

%% downsample to 200 Hz defualt sampling rate of Pan-tompkins for those Fs
%% which are devidable 
out =~rem(fs,200)*fs/200 ;
if out ~= 0
ecg = downsample(ecg,out);
fs = 200; %new fs
end

%% Plot differently based on filtering settings
% if fs == 200
%  figure,  ax(1)=subplot(321);plot(ecg);axis tight;title('Raw ECG Signal');
% else
%  figure,  ax(1)=subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw ECG Signal');
% end
    
%% Noise cancelation(Filtering) % Filters (Filter in between 5-15 Hz)
if fs == 200
%% Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
h_l = filter(b,a,[1 zeros(1,12)]); 
ecg_l = conv (ecg ,h_l);
ecg_l = ecg_l/ max( abs(ecg_l));
delay = 6; %based on the paper
%ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');
%% High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];
h_h = filter(b,a,[1 zeros(1,32)]); 
ecg_h = conv (ecg_l ,h_h);
ecg_h = ecg_h/ max( abs(ecg_h));
delay = delay + 16; % 16 samples for highpass filtering
%ax(3)=subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');
else
%% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
f1=5; %cuttoff low frequency to get rid of baseline wander
f2=15; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg_h = filtfilt(a,b,ecg);
ecg_h = ecg_h/ max( abs(ecg_h));
%ax(3)=subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
end
%% derivative filter H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
h_d = [-1 -2 0 2 1]*(1/8);%1/8*fs
ecg_d = conv (ecg_h ,h_d);
ecg_d = ecg_d/max(ecg_d);
delay = delay + 2; % delay of derivative filter 2 samples
%ax(4)=subplot(324);plot(ecg_d);axis tight;title('Filtered with the derivative filter');
%% Squaring nonlinearly enhance the dominant peaks
ecg_s = ecg_d.^2;
%ax(5)=subplot(325);plot(ecg_s);axis tight;title('Squared');




%% Moving average Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + 15;
%ax(6)=subplot(326);plot(ecg_m);axis tight;title('Averaged with 30 samples length,Black noise adaptive threshold,Green QRS adaptive threshold');
%axis tight;


%% Fiducial Mark 
% Note : a minimum distance of 40 samples is considered between each R wave
% since in physiological point of view no RR wave can occur in less than
% 200 msec distance
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',0.2*fs);


%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(ecg_m(1:2*fs))*1/4; % 0.25 of the max amplitude 
THR_NOISE = mean(ecg_m(1:2*fs))*1/2; % 0.5 of the mean signal is considered to be noise

%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
THR_SIG1 = max(ecg_h(1:2*fs))*1/4; % 0.25 of the max amplitude 
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2; %
%% Thresholding and online desicion rule

for i = 1 : length(pks)
    
    
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)      % Fixed from 0.100 to 0.150
          [y_i x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
       else
          if i == 1
            [y_i x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
          elseif locs(i)>= length(ecg_h)
            [y_i x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
          end
        
     end
    
    
  % update the heart_rate  
    if length(qrs_c) >= 9 
          
          diffRR = diff(qrs_i(end-8:end)); %calculate RR interval
          mean_RR = mean(diffRR); % calculate the mean of 8 previous R waves interval
          temp = diffRR(diffRR >= 0.92*mean_RR); %find the most regular beats
          temp = temp(temp <= 1.16*mean_RR);
          if isempty(temp)
              % lower down thresholds to detect better in MVI
               THR_SIG = 0.5*(THR_SIG);
               THR_NOISE = 0.5*(THR_NOISE);  
              % lower down thresholds to detect better in Bandpass filtered 
               THR_SIG1 = 0.5*(THR_SIG1);
               THR_NOISE1 = 0.5*(THR_NOISE1); 
          end
          selected_RR =[selected_RR temp];
          % keep a track of the most regular RR intervals 
           if  length(selected_RR) >= 8
              m_selected_RR = mean(selected_RR);
              
              if length(selected_RR)-8 > 0
                selected_RR(1:length(selected_RR)-8)=[]; 
              end
           end   
    end
    
    % find noise and QRS peaks
    if pks(i) >= THR_SIG
        
                 % if a QRS candidate occurs within 360ms of the previous QRS
                 % ,the algorithm determines if its T wave or QRS
                 if length(qrs_c) >= 3
                      if (locs(i)-qrs_i(end)) <= round(0.3600*fs)
                        Slope1 = mean(diff(ecg_m(locs(i)-round(0.050*fs):locs(i)))); %mean slope of the waveform at that position
                        Slope2 = mean(diff(ecg_m(qrs_i(end)-round(0.050*fs):qrs_i(end)))); %mean slope of previous R wave
                             if abs(Slope1) <= abs(0.5*(Slope2)) || (locs(i)-qrs_i(end)) <= round(0.4*test_m)  % slope less then 0.5 of previous R
                                 nois_c = [nois_c pks(i)];
                                 nois_i = [nois_i locs(i)];
                                 skip = 1; % T wave identification
                             else
                                 skip = 0;
                             end
            
                      end
                 end
        
        if skip == 0  % skip is 1 when a T wave is detected       
        qrs_c = [qrs_c pks(i)];
        qrs_i = [qrs_i locs(i)];
        
        % bandpass filter check threshold
         if y_i >= THR_SIG1
                        if ser_back 
                           qrs_i_raw = [qrs_i_raw x_i];  % save index of bandpass 
                        else
                           qrs_i_raw = [qrs_i_raw locs(i)-round(0.150*fs)+ (x_i - 1)];% save index of bandpass 
                        end
                           qrs_amp_raw =[qrs_amp_raw y_i];% save amplitude of bandpass 
          SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;% adjust threshold for bandpass filtered sig
         end
         
        % adjust Signal level
        SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;
        end
        
        
    elseif THR_NOISE <= pks(i) && pks(i)<THR_SIG
        
       % calculate the mean of the last 8 R waves to make sure that QRS is not
       % missing(If no R detected , trigger a search back) 1.66*mean
       
       if m_selected_RR
           test_m = m_selected_RR; %if the regular RR availabe use it
           
       elseif mean_RR && m_selected_RR == 0
           test_m = mean_RR;
           
       else
           test_m = 0;
       end
        
    if test_m
          if (locs(i) - qrs_i(end)) >= round(1.66*test_m)% it shows a QRS is missed 
              pks(i) = max(ecg_m(qrs_i(end)+ round(0.200*fs):locs(i))); % search back and locate the max in this interval
             if pks(i) >= THR_NOISE
               qrs_c = [qrs_c pks(i)];
               qrs_i = [qrs_i locs(i)];
               
               % take care of bandpass signal threshold
               if y_i >= THR_NOISE1 
                        if ser_back 
                           qrs_i_raw = [qrs_i_raw x_i];  % save index of bandpass 
                        else
                           qrs_i_raw = [qrs_i_raw locs(i)-round(0.150*fs)+ (x_i - 1)];% save index of bandpass ...  %% Fixed from 0.*fs to 0.150*fs
                        end
                           qrs_amp_raw =[qrs_amp_raw y_i]; %save amplitude of bandpass 
                           SIG_LEV1 = 0.25*y_i + 0.75*SIG_LEV1; %when found with the second thres 
               end
               
               not_nois = 1;
               SIG_LEV = 0.25*pks(i) + 0.75*SIG_LEV ;  %when found with the second threshold             
             end 
              
          else
              not_nois = 0;
              
          end
    end
       
    
      if not_nois == 0 %it is not noise when no_nois =1 
        nois_c = [nois_c pks(i)];
        nois_i = [nois_i locs(i)];
        
        
        if THR_NOISE1 <= y_i %&&  y_i < THR_SIG1
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        end
        
         %adjust Noise level
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;  
        
        
      end
        
           
    end
    
    % adjust the threshold with SNR
    if NOISE_LEV ~= 0 && SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.2*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    % adjust the threshold with SNR for bandpassed signal
    if NOISE_LEV1 ~= 0 && SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.2*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    
    
    
    
 skip = 0; %reset parameters
 not_nois = 0; %reset parameters
 ser_back = 0;  %reset bandpass param   
end
%hold on,scatter(qrs_i,qrs_c,'r');
%hold on,plot(ones(1,length(ecg_m))*THR_NOISE,'--k');
%hold on,plot(ones(1,length(ecg_m))*THR_SIG,'-.r');
%if ax(:)
%linkaxes(ax,'x');
%zoom on;
%end



%% overlay on the signals
%figure,az(1)=subplot(311);plot(ecg_h);title('QRS on Filtered Signal');axis tight;
%hold on,scatter(qrs_i_raw,qrs_amp_raw,'r');
%az(2)=subplot(312);plot(ecg_m);title('QRS on MVI signal and Noise level(green), and Threshold(black)');axis tight;
%hold on,scatter(qrs_i,qrs_c,'r');
%hold on,plot(ones(1,length(ecg_m))*THR_NOISE,'LineWidth',2,'Linestyle','--','color','g');
%hold on,plot(ones(1,length(ecg_m))*THR_SIG,'LineWidth',2,'Linestyle','-.','color','k');
%az(3)=subplot(313);plot(ecg-mean(ecg));title('Pulse train of the found QRS on ECG signal');axis tight;
%line(repmat(qrs_i_raw,[2 1]),repmat([min(ecg-mean(ecg))/2; max(ecg-mean(ecg))/2],size(qrs_i_raw)),'LineWidth',2,'LineStyle','-.','Color','r'),
%linkaxes(az,'x');
%zoom on;

end
 









