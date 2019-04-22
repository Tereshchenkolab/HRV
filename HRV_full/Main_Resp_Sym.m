%Respiration

clear all
close all
clc

warning('off');


% load mat files
[fname,pname] = uigetfile('*','Verzeichnis wählen');
tachofiles = dir(pname);
anz=length(tachofiles);

pname_saving='C:\Users\kabir\Desktop\New folder\';         %Copy Folder


for fileCT=3:1:anz %1
    
  close all;
  
  try
           
      filename=tachofiles(fileCT).name;
      file=strcat(pname,filename);
      [pathstr, full_name, ext] = fileparts(file);
    
      [name, remain] = strtok(full_name, '_');
    
    
      disp(strcat(pname,name));
    
    
      FidPts=[];
      FidPts=load(strcat(pname,full_name,'.mat'));
      
      fs=[]; xi=[]; yi=[]; zi=[];
      r_points_x=[]; r_points_y=[]; r_points_z=[];
      q_points_x=[]; q_points_y=[]; q_points_z=[];
      s_points_x=[]; s_points_y=[]; s_points_z=[];
      
      fs=FidPts.sampling_rate;
      xi=FidPts.xi;
      yi=FidPts.yi;
      zi=FidPts.zi;
      
      r_points_x=FidPts.r_points_x;
      r_points_y=FidPts.r_points_y;
      r_points_z=FidPts.r_points_z;
      
      q_points_x=FidPts.q_points_x;
      q_points_y=FidPts.q_points_y;
      q_points_z=FidPts.q_points_z;
      
      s_points_x=FidPts.s_points_x;
      s_points_y=FidPts.s_points_y;
      s_points_z=FidPts.s_points_z;
      
      LD=1; % X lead
%       samp_3min=3*60*fs; % samples in 3 minutes
      
      samp_10s=10*fs; % samples in 10 sec
      samp_3min=samp_10s;
      
      R_samp3_3min=[];
      Q_samp3_3min=[];
      S_samp3_3min=[];
      XYZ_samp3_3min=[];
      
      if LD==1
          
%           %3 minute
%           st_samp=[]; ed_samp=[];
%           st_samp=q_points_x(2,1)-5; %start sample; start from 2nd Qpoint minus 5 samples
%           ed_samp=min((st_samp+samp_3min),length(xi));
%           locs_Spts=[];
%           locs_Spts=find(s_points_x(:,1)>=st_samp & s_points_x(:,1)<=ed_samp);
%           S_samp3_3min=s_points_x(locs_Spts,1);          
%           locs_Qpts=[];
%           locs_Qpts=find(q_points_x(:,1)>=st_samp & q_points_x(:,1)<S_samp3_3min(end));
%           Q_samp3_3min=q_points_x(locs_Qpts,1);          
%           locs_Rpts=[];
%           locs_Rpts=find(r_points_x(:,1)>=st_samp & r_points_x(:,1)<S_samp3_3min(end));
%           R_samp3_3min=r_points_x(locs_Rpts,1);
          
          %10s - 2nd part
          st_samp=[]; ed_samp=[];
          ed_samp=s_points_x(end-1,1)+5; %end sample; end at 2nd-last Rpeak plus 5 samples
          st_samp=ed_samp-samp_10s;
          locs_Qpts=[];
          locs_Qpts=find(q_points_x(:,1)>=st_samp & q_points_x(:,1)<=ed_samp);
          Q_samp3_3min=q_points_x(locs_Qpts,1);
          locs_Rpts=[];
          locs_Rpts=find(r_points_x(:,1)>=Q_samp3_3min(1) & r_points_x(:,1)<=ed_samp);
          R_samp3_3min=r_points_x(locs_Rpts,1);
          locs_Spts=[];
          locs_Spts=find(s_points_x(:,1)>=Q_samp3_3min(1) & s_points_x(:,1)<=ed_samp);
          S_samp3_3min=s_points_x(locs_Spts,1);  
          
          XYZ_samp3_3min=xi(st_samp:ed_samp);
          S_samp3_3min=S_samp3_3min-st_samp+1;
          Q_samp3_3min=Q_samp3_3min-st_samp+1;
          R_samp3_3min=R_samp3_3min-st_samp+1;
          
      elseif LD==2
%           
%           %3 minute
%           st_samp=[]; ed_samp=[];
%           st_samp=q_points_y(2,1)-5; %start sample; start from 2nd Rpeak minus 5 samples
%           ed_samp=min((st_samp+samp_3min),length(yi));
%           locs_Spts=[];
%           locs_Spts=find(s_points_y(:,1)>=st_samp & s_points_y(:,1)<=ed_samp);
%           S_samp3_3min=s_points_y(locs_Spts,1);         
%           locs_Qpts=[];
%           locs_Qpts=find(q_points_y(:,1)>=st_samp & q_points_y(:,1)<S_samp3_3min(end));
%           Q_samp3_3min=q_points_y(locs_Qpts,1);         
%           locs_Rpts=[];
%           locs_Rpts=find(r_points_y(:,1)>=st_samp & r_points_y(:,1)<S_samp3_3min(end));
%           R_samp3_3min=r_points_y(locs_Rpts,1);
          
          %10s - 2nd part
          st_samp=[]; ed_samp=[];
          ed_samp=s_points_y(end-1,1)+5; %end sample; end at 2nd-last Rpeak plus 5 samples
          st_samp=ed_samp-samp_10s;
          locs_Qpts=[];
          locs_Qpts=find(q_points_y(:,1)>=st_samp & q_points_y(:,1)<=ed_samp);
          Q_samp3_3min=q_points_y(locs_Qpts,1);
          locs_Rpts=[];
          locs_Rpts=find(r_points_y(:,1)>=Q_samp3_3min(1) & r_points_y(:,1)<=ed_samp);
          R_samp3_3min=r_points_y(locs_Rpts,1);
          locs_Spts=[];
          locs_Spts=find(s_points_y(:,1)>=Q_samp3_3min(1) & s_points_y(:,1)<=ed_samp);
          S_samp3_3min=s_points_y(locs_Spts,1);
          
          XYZ_samp3_3min=yi(st_samp:ed_samp);
          S_samp3_3min=S_samp3_3min-st_samp+1;
          Q_samp3_3min=Q_samp3_3min-st_samp+1;
          R_samp3_3min=R_samp3_3min-st_samp+1;
          
      elseif LD==3
          
%           %3 minute
%           st_samp=[]; ed_samp=[];
%           st_samp=r_points_z(2,1)-5; %start sample; start from 2nd Rpeak minus 5 samples
%           ed_samp=min((st_samp+samp_3min),length(zi));
%           locs_Spts=[];
%           locs_Spts=find(s_points_z(:,1)>=st_samp & s_points_z(:,1)<=ed_samp);
%           S_samp3_3min=s_points_z(locs_Spts,1);         
%           locs_Qpts=[];
%           locs_Qpts=find(q_points_z(:,1)>=st_samp & q_points_z(:,1)<S_samp3_3min(end));
%           Q_samp3_3min=q_points_z(locs_Qpts,1);         
%           locs_Rpts=[];
%           locs_Rpts=find(r_points_z(:,1)>=st_samp & r_points_z(:,1)<S_samp3_3min(end));
%           R_samp3_3min=r_points_z(locs_Rpts,1);
          
          %10s - 2nd part
          st_samp=[]; ed_samp=[];
          ed_samp=s_points_z(end-1,1)+5; %end sample; end at 2nd-last Rpeak plus 5 samples
          st_samp=ed_samp-samp_10s;
          locs_Qpts=[];
          locs_Qpts=find(q_points_z(:,1)>=st_samp & q_points_z(:,1)<=ed_samp);
          Q_samp3_3min=q_points_z(locs_Qpts,1);
          locs_Rpts=[];
          locs_Rpts=find(r_points_z(:,1)>=Q_samp3_3min(1) & r_points_z(:,1)<=ed_samp);
          R_samp3_3min=r_points_z(locs_Rpts,1);
          locs_Spts=[];
          locs_Spts=find(s_points_z(:,1)>=Q_samp3_3min(1) & s_points_z(:,1)<=ed_samp);
          S_samp3_3min=s_points_z(locs_Spts,1);
          
          XYZ_samp3_3min=zi(st_samp:ed_samp);
          S_samp3_3min=S_samp3_3min-st_samp+1;
          Q_samp3_3min=Q_samp3_3min-st_samp+1;
          R_samp3_3min=R_samp3_3min-st_samp+1;
          
      end
      
      
      Rms_samp3_3min=[]; %R in ms
      Rms_samp3_3min=(R_samp3_3min/fs)*1000;
      
      
      
      %% Respiration Analysis
        
           %8ms
            samp8=(8/1000)*fs;

            % figure; plot(xi);hold on;plot(r_points_x(:,1),r_points_x(:,2),'.r');plot(q_points_x(:,1),q_points_x(:,2),'.g');plot(s_points_x(:,1),s_points_x(:,2),'.k');hold off


            % 4Hz interpolation
            fs_interp=4;


            %Using VCG-X
            resp_points_US_X=[];
            resp_points_DS_X=[];
            resp_points_RA_X=[];

            [resp_points_US_X,resp_points_DS_X,resp_points_RA_X]=resp_points(XYZ_samp3_3min,Q_samp3_3min,R_samp3_3min,S_samp3_3min,fs,samp8);

            time_X=[];
            resp_interp_US_X=[]; resp_interp_DS_X=[]; resp_interp_RA_X=[];
        %     tt=1/fs_interp:1/fs_interp:resp_points_US(end,1)/fs;
            time_X=0:1/fs_interp:length(XYZ_samp3_3min)/fs;
            resp_interp_US_X=spline(resp_points_US_X(:,1)/fs,resp_points_US_X(:,2),time_X);
            resp_interp_DS_X=spline(resp_points_DS_X(:,1)/fs,resp_points_DS_X(:,2),time_X);
            resp_interp_RA_X=spline(resp_points_RA_X(:,1)/fs,resp_points_RA_X(:,2),time_X);


            %Original signals
            T_X=[]; %time
            T_X=1/fs:1/fs:length(XYZ_samp3_3min)/fs;
            
%             RespFig=figure('Visible','off');
%             plot(T_X,XYZ_samp3_3min);hold on;plot(time_X,resp_interp_US_X.*10,'r');plot(time_X,resp_interp_DS_X.*10,'g');plot(time_X,resp_interp_RA_X.*10,'k');hold off
            

            %Same sampling points as original signal 
            tempT=[]; tempX=[];
            tempT=1/fs:1/fs:length(XYZ_samp3_3min)/fs;
            tempX=interp1(time_X,resp_interp_US_X,tempT);
            
            tempX(isnan(tempX))=0;
            
            RespSig=[];
            RespSig=tempX;
            
%             %Polynomial Fitting
%             Pdeg=20;
%             P_L=[];
%             P_L=polyfit(1:length(tempX),tempX,Pdeg);
%         
%             RespSig=[];
%             RespSig=polyval(P_L,1:length(tempX));
            
            
            mag_X_R=[]; 
            mag_X_R=RespSig(R_samp3_3min);
            
            
            RespFig=figure('Visible','off');
            plot(T_X,XYZ_samp3_3min);hold on;plot(tempT,RespSig.*10,'r');hold off
            
            
            %Hilbert Transform

            ht_X=[]; ph_X=[]; ph_X_R=[]; 
            ht_X=hilbert(RespSig);
            ph_X=angle(ht_X);
            ph_X_R=ph_X(R_samp3_3min);
            

%             [maxValue_X,indexMax_X] = max(abs(fft(RespSig-mean(RespSig))));
%             resp_frequency_X = indexMax_X * fs / length(RespSig);
% 
%             Ax_fil=[];
%             Wn=resp_frequency_X/(fs/2); % cutt off based on fs
%             N = 3; % order of 3 less processing
%             [b,a] = butter(N,Wn,'high'); %bandpass filtering
%             Ax_fil=filtfilt(b,a,xi);


            
      
      %% Analysis
      
            RRms_All=[]; BeatCT=0;          
            RRms_All=diff(Rms_samp3_3min);
            BeatCT=length(Rms_samp3_3min);
            
            ph_X_RRms=[];
            ph_X_RRms=ph_X_R(2:end); % To match respiratory phase with RR
            

      %% Symbolic Analysis

            %Aalso includes analysis for Alternans excluding respiratory phase transitions
            
            Mat_mk=NaN(5,6);



                % Threshold 
                threshold=[0 1 2 3 4]; %in samples
                threshold=(threshold/fs)*1000; %in ms

                for sdCT=1:length(threshold)

                    SDmk_R0=0; SDmk_R1=0; SDmk_LR2=0; SDmk_ULR2=0; SDmk_ALT=0; SDmk_ALT_exResTr=0; SDmk_allSym=0; 


                        
                        RRms_All_Th=[];
                        if sdCT==1
                            RRms_All_Th=RRms_All;
                        else
                            RRms_All_Th=round(RRms_All/threshold(sdCT))*threshold(sdCT);
                        end
                        
                        [SDmk,allSym]=calc_CarResSym(RRms_All_Th,ph_X_RRms);
                        

                        SDmk_R0=SDmk_R0+SDmk.R0;
                        SDmk_R1=SDmk_R1+SDmk.R1;
                        SDmk_LR2=SDmk_LR2+SDmk.LR2;
                        SDmk_ULR2=SDmk_ULR2+SDmk.ULR2;
                        SDmk_ALT=SDmk_ALT+SDmk.ALT;
                        SDmk_ALT_exResTr=SDmk_ALT_exResTr+SDmk.ALT_exResTr;
                        SDmk_allSym=SDmk_allSym+allSym;


                    Mat_mk(sdCT,1)=SDmk_R0/SDmk_allSym*100;
                    Mat_mk(sdCT,2)=SDmk_R1/SDmk_allSym*100;
                    Mat_mk(sdCT,3)=SDmk_LR2/SDmk_allSym*100;
                    Mat_mk(sdCT,4)=SDmk_ULR2/SDmk_allSym*100;
                    Mat_mk(sdCT,5)=SDmk_ALT/SDmk_allSym*100;
                    Mat_mk(sdCT,6)=SDmk_ALT_exResTr/SDmk_allSym*100;

                end
      
            
            

      %%  Save            

           
%             excelName='HRV_ALTSymbolic_Results_3min.xls';
            excelName='HRV_ALTSymbolic_Results_10s_2.xls';
            
                
            if ~exist(strcat(pname_saving,excelName))
                fid1 = fopen(strcat(pname_saving,excelName),'a');
                fprintf(fid1,'File ID \t');
                fprintf(fid1,'Number of beats analysed \t');

                fprintf(fid1,' R0 Th0 \t');
                fprintf(fid1,' R1 Th0 \t');
                fprintf(fid1,' LR2 Th0 \t');
                fprintf(fid1,' ULR2 Th0 \t');
                fprintf(fid1,' ALT Th0 \t');
                fprintf(fid1,' ALT_exResTr Th0 \t');
                %
                fprintf(fid1,' R0 Th1 \t');
                fprintf(fid1,' R1 Th1 \t');
                fprintf(fid1,' LR2 Th1 \t');
                fprintf(fid1,' ULR2 Th1 \t');
                fprintf(fid1,' ALT Th1 \t');
                fprintf(fid1,' ALT_exResTr Th1 \t');
                %
                fprintf(fid1,' R0 Th2 \t');
                fprintf(fid1,' R1 Th2 \t');
                fprintf(fid1,' LR2 Th2 \t');
                fprintf(fid1,' ULR2 Th2 \t');
                fprintf(fid1,' ALT Th2 \t');
                fprintf(fid1,' ALT_exResTr Th2 \t');
                %
                fprintf(fid1,' R0 Th3 \t');
                fprintf(fid1,' R1 Th3 \t');
                fprintf(fid1,' LR2 Th3 \t');
                fprintf(fid1,' ULR2 Th3 \t');
                fprintf(fid1,' ALT Th3 \t');
                fprintf(fid1,' ALT_exResTr Th3 \t');
                %
                fprintf(fid1,' R0 Th4 \t');
                fprintf(fid1,' R1 Th4 \t');
                fprintf(fid1,' LR2 Th4 \t');
                fprintf(fid1,' ULR2 Th4 \t');
                fprintf(fid1,' ALT Th4 \t');
                fprintf(fid1,' ALT_exResTr Th4 \t');

                fprintf(fid1,'\n \n');
                fclose(fid1);

            end



            fid1 = fopen(strcat(pname_saving,excelName),'a');
            fprintf(fid1,'%s \t',name);
            fprintf(fid1,'%d \t ',BeatCT); 

            for k1=1:length(Mat_mk(:,1))
                for k2=1:length(Mat_mk(1,:))
                    if ~isnan(Mat_mk(k1,k2))
                        fprintf(fid1,'%2.1f \t ',Mat_mk(k1,k2));
                    else
                        fprintf(fid1,'\t ');
                    end
                end
            end

            fprintf(fid1,'\n');
            fclose(fid1);
            
          
                
                
            % Save Figures

%             figName='_3min.jpg';
            figName='_10s_2.jpg';
            
%             Images_folder = [pname_saving 'All Figures Respiration' '\'];
            Images_folder = [pname_saving 'All Figures Respiration 10s_2' '\'];
            if (exist(Images_folder) == 0)
                mkdir (Images_folder);
            end
           
            filepath = [Images_folder strcat(name,figName)];
            saveas(RespFig,filepath,'jpg');



%             Mat_folder = [pname_saving 'Mat Files Respiration' '\'];
            Mat_folder = [pname_saving 'Mat Files Respiration 10s_2' '\'];
            if (exist(Mat_folder) == 0)
                mkdir (Mat_folder);
            end

            save(strcat(Mat_folder,'\',name,'.mat'), 'R_samp3_3min','Q_samp3_3min','S_samp3_3min','XYZ_samp3_3min','LD','samp_3min','fs',...
                'time_X','resp_interp_US_X','resp_interp_DS_X','resp_interp_RA_X','tempT','tempX','RespSig','ht_X','ph_X','ph_X_R');   

            
    
  catch
      
%       disp('error');
      
      fid1 = fopen(strcat(pname_saving,'Removed.xls'),'a');
      fprintf(fid1,'%s \t',name);
      fprintf(fid1,'\n');
      fclose(fid1);
        
  end


end