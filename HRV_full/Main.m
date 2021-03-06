%Analysis for HRV

clear all
close all
clc

warning('off');

fs = 200; %sampling rate of data


folder_root = ''; %%file path to data folder

folder_dir = dir(strcat(folder_root,'13*'));

for fold_id=1:length(folder_dir)

    %File path name for saving
    pname_saving=strcat(folder_root,folder_dir(fold_id).name,'/');

    Reviewed_folder = [pname_saving 'Poincare Plots' '/'];
    if (exist(Reviewed_folder) == 0)
      mkdir (Reviewed_folder);
    end

    tachofiles=dir(strcat(pname_saving,'1*.mat'));

    excelName=strcat(folder_root,'Excel_results/',folder_dir(fold_id).name,'_HRV_3min.xls');
    
    fid1 = fopen(strcat(excelName),'a');
    fprintf(fid1,'File ID \t');
    fprintf(fid1,'Number of Hour (hr) \t');
    fprintf(fid1,'Time of 3 min (min) \t');
    fprintf(fid1,'Total of beats in the 3 mins \t');

    fprintf(fid1,' Mean RR (ms) \t');
    fprintf(fid1,' RMSSD (ms) \t');

    fprintf(fid1,' LF \t');
    fprintf(fid1,' HF \t');
    fprintf(fid1,' LF/HF \t');

    fprintf(fid1,' Mean relative RR (ms) \t');
    fprintf(fid1,' SD relative RR (ms) \t');
    fprintf(fid1,' relative HRV (ms) \t');

    fprintf(fid1,' Sample Entropy \t');
    fprintf(fid1,' Renyi Entropy \t');

    fprintf(fid1,' Poincare plot SD1 \t');
    fprintf(fid1,' Poincare plot SD2 \t');
    fprintf(fid1,' Poincare plot SD1/SD2 \t');

    fprintf(fid1,'\n \n');
    fclose(fid1);
  

    anz=length(tachofiles);


    for fileCT=1:anz 
 
        filename=tachofiles(fileCT).name;
        file=strcat(pname_saving,filename);
    
        C_s=strsplit(filename,'_');
        file_id=strcat(C_s{1},'_',C_s{2});

        hour=C_s{length(C_s)-1};
        
        FidPts=[];
        r_points_x=[];
        time_3min=[];
        FidPts=[];
        I_d=[];
        M_d=[];

        FidPts=load(file);
    
        M_d=FidPts.good_segment;

        if isempty(M_d)==0 


            r_points_x=(FidPts.good_Rpeaks/fs)*1000; %converts good_Rpeaks from samples (calculated by r-peak main code) to ms 
 
            RRinterval=diff(r_points_x); %calculates RR intervals from good_Rpeaks in ms 
 
 
            %%samp_3min=3*60*fs; % number of samples in 3 minutes
      
            R_samp3_3min=[];
          
            R_samp3_3min=r_points_x;

            Rms_samp3_3min=R_samp3_3min;
	
            Total_beats=[];
            Total_beats=length(Rms_samp3_3min);      
    
            time_3min=(FidPts.time_3min)*3
   
  %% Analysis

      %calculate RR intervals  
            RRms_All=[]; BeatCT=0;
            segCT=3;
            RRms_All=diff(Rms_samp3_3min); %RR intervals in ms
            BeatCT=length(Rms_samp3_3min);


            % Basic parameters
            meanRR=[]; sdRR=[]; RMSSD=[]; pNN50=[];
            meanRR=mean(RRms_All(~isnan(RRms_All)));
            sdRR=std(RRms_All(~isnan(RRms_All)));
            RMSSD=sqrt((sum((diff(RRms_All)).^2))/(length(RRms_All)-1));
            pNN50=(sum(abs(diff(RRms_All))>50)/(length(RRms_All)-1))*100;
            
            % Spectral analysis
            RRs_All=[];
            RRs_All=RRms_All./1000;
            resample_fs = 500; % 4; 250; 1000
            [LF,HF,LFHFratio] = spectral_analysis_HRV(RRs_All,resample_fs);

            % Relative RR
            rel_RR=[]; rel_RR_HRV=[];
            [rel_RR,rel_RR_HRV]=relativeRR(RRms_All);

            %Sample Entropy
            sampentr = [];
            dim = 2;                  % embedded dimension
            r = 0.2 * std(RRms_All);  % tolerance (typically 0.2 * std)
            tau = 1;                  % delay time for downsampling
            sampentr = SampEn( dim, r, RRms_All(~isnan(RRms_All)), tau );

            %Shannon and Renyi Entropy
            shannonEntr4=[]; shannonEntr3=[]; renyiEntr4=[];
            byteLevel4=4; %4ms window
            byteLevel3=3; %3ms window
            % alpha = 3 or 4 for Renyi entropy improved accuracy 
            % in separating congestive heart failure from normal sinus rhythm.
            % Conforth et al., CinC 2016
            alpha=4;
            [shannonEntr4, renyiEntr4] = shanRenEntropy(RRms_All(~isnan(RRms_All)),byteLevel4,alpha);
            [shannonEntr3, ~] = shanRenEntropy(RRms_All(~isnan(RRms_All)),byteLevel3,alpha);

            %Poincare Plot
            PP_SD1=[]; PP_SD2=[]; PP_S1S2=[];
            PP_SD1=sqrt(var(diff(RRms_All)/sqrt(2)));
            for i=1:length(RRms_All)-1
                s_dRRms(i)=RRms_All(i)+RRms_All(i+1);
            end
            PP_SD2=sqrt((2*(sdRR^2))-((PP_SD1)^2));
            PP_S1S2=PP_SD1/PP_SD2;
            
            %Save poincare figures
            fig=figure('visible','off');
            plot(RRms_All(1:end-1),RRms_All(2:end),'ko');
            saveas(gcf,strcat(Reviewed_folder,'\',folder_dir(fold_id).name,'_HOUR_',hour,'.jpg'),'jpg');

            %%  Save            

            excelName=strcat(folder_root,'Excel_results/',folder_dir(fold_id).name,'_HRV_3min.xls');
            fid1 = fopen(strcat(excelName),'a');

	        fprintf(fid1,'%s \t',file_id);
        	fprintf(fid1,'%s \t',hour);
        	fprintf(fid1,'%d \t',time_3min);
            fprintf(fid1,'%d \t',Total_beats);

            fprintf(fid1,'%2.3f \t ',meanRR);
            fprintf(fid1,'%2.3f \t ',RMSSD);

            fprintf(fid1,'%2.3f \t ',LF);
            fprintf(fid1,'%2.3f \t ',HF);
            fprintf(fid1,'%2.3f \t ',LFHFratio);

            if ~isempty(rel_RR)
                fprintf(fid1,'%2.3f \t ',mean(rel_RR(~isnan(rel_RR))));
                fprintf(fid1,'%2.3f \t ',std(rel_RR(~isnan(rel_RR))));
            else
                fprintf(fid1,'\t \t');
            end
            if ~isempty(rel_RR_HRV) && ~isnan(rel_RR_HRV)
                fprintf(fid1,'%2.3f \t ',rel_RR_HRV);
            else
                fprintf(fid1,'\t ');
            end

            if ~isempty(sampentr) && ~isnan(sampentr) && ~isinf(sampentr)
                fprintf(fid1,'%2.3f \t ',sampentr);
            else
                fprintf(fid1,'\t ');
            end

            if ~isempty(renyiEntr4) && ~isnan(renyiEntr4)
                fprintf(fid1,'%2.3f \t ',renyiEntr4);
            else
                fprintf(fid1,'\t ');
            end

            fprintf(fid1,'%2.3f \t ',PP_SD1);
            fprintf(fid1,'%2.3f \t ',PP_SD2);
            fprintf(fid1,'%2.3f \t ',PP_S1S2);
            
            fprintf(fid1,'\n \n');
            fclose(fid1);

         
        else

%         excelName='HRV_Results_3min.xls';

%         if ~exist(strcat(pname_saving,excelName))

            excelName=strcat(folder_root,'Excel_results/',folder_dir(fold_id).name,'_HRV_3min.xls');
            fid1 = fopen(strcat(excelName),'a');
%            fid1 = fopen(strcat(pname_saving,excelName),'a');

            fprintf(fid1,'%s \t',file_id);
            fprintf(fid1,'%s \t',hour);
            fprintf(fid1,'\n');
            fclose(fid1);

 %	     end

        
        end

    end

    close all;

end
