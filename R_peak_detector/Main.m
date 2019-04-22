%This is the main code to detect the R-peaks in a 3 minute sinus epoch in each hour of data. 
%First it will run a modified pan-tompkins algorithm to
%detect all r-peaks, and then it will look at the first 3-minute window of
%signal to detect whether the duration of each RR interval is within 15% of
%the previous interval. If it is, the algorithm selects this 3 minute epoch
%and labels the good_Rpeaks for use in the next code for HRV analysis. If
%not, the window slides to the next 3-minute segment after the bad interval
%and recalculates duration. If an entire hour is scanned
%and no good segment is found, the algorithm moves to the next r-peak
%detection algorithm and tries again to ensure proper r-peak calculation.
clear all
close all

%Read root
fs=200; %sampling rate of data

folder_root='X:\OHSU Shared\Restricted\SOM\MED\CARDIO\Tereshchenko\PACE Study\HRV Study (Nichole)\Code\Code for GitHub\Sample Data\' %add file path for data

excel_root='X:\OHSU Shared\Restricted\SOM\MED\CARDIO\Tereshchenko\PACE Study\HRV Study (Nichole)\Code\Code for GitHub\Sample Data\Excel_results\'; %add file path for excel files to be saved after calculation

tachofolder=dir(strcat(folder_root,'1*Z*'));

anz=length(tachofolder);

% start parallel loop
parfor folder_dat=1:anz

	file_root=strcat(folder_root,tachofolder(folder_dat).name,'\');

	tachofiles=dir(strcat(file_root,'1*.mat'));	

	anz2=length(tachofiles);



    % create an excel file for each case
	excelName=strcat(excel_root,tachofolder(folder_dat).name,'.xls');

        fid1 = fopen(strcat(excelName),'a');
        fprintf(fid1,'File ID \t');
        fprintf(fid1,'Number of Hour (hr) \t');
        fprintf(fid1,'Time of 3 min (min) \t');
        fprintf(fid1,'Total of beats in the 3 mins \t');
        fprintf(fid1,' Mean RR (ms) \t');
        fprintf(fid1,' SD RR (ms) \t');
        fprintf(fid1,' RMSSD (ms) \t');
        fprintf(fid1,' Algorithm Number \t');
        fprintf(fid1,' Number over 15 Percent \t'); %These "number over" percents are for review to identify potential problems in detecting sinus rhythm
        fprintf(fid1,' Number over 12 Percent \t'); %They dont do anything to the calculation
        fprintf(fid1,' Number over 10 Percent \t');
        fprintf(fid1,' Number over 8 Percent \t');
        fprintf(fid1,' Number over 6 Percent \t');
        fprintf(fid1,' Number over 4 Percent \t');
        fprintf(fid1,' Number over 2 Percent \t');

        fprintf(fid1,'\n \n');
        fclose(fid1);

    % read each file
	for file_dat=1:anz2




		file_name=tachofiles(file_dat).name

	        C_s=strsplit(file_name,'_');
        	file_id=strcat(C_s{1},'_',C_s{2});
	        hour=C_s{end-1};

	% load the file and check if the data exist
		matfile=load(strcat(file_root,file_name));

		if isfield(matfile,'bw_data')
            sequence_length = length(matfile.bw_data);
			

			peakloc=[];
			data=[];
			pos_d=[];
                 	time_3min=[];
                        Total_beats=[];

                        meanRR=[];
                        sdRR=[];
                        RMSSD=[];

            % read each 3 minutes segment 
            algo_num = 1;
            interval_t = 36000;
            good_percent = 0.15; % threshold for determining sinus
            check_percent = 0.15;
            check_percent2 = 0.12;
            check_percent3 = 0.10;
            check_percent4 = 0.08;
            check_percent5 = 0.06;
            check_percent6 = 0.04;
            check_percent7 = 0.02;
            
            while true 
              pos_start = 1;
              while pos_start+interval_t < sequence_length

                start_index=pos_start;
                end_index=pos_start + interval_t - 1;
                disp('start_index:'), disp(start_index)
                disp('end_index:'), disp(end_index)

                data=matfile.bw_data(start_index:end_index);

                if algo_num == 1
                  peakloc=rpeak_pan_tompkins(data,fs);
                elseif algo_num == 2
                  peakloc=rpeak_parabolic(data,fs);
                else
                  peakloc=rpeak_pca(data,fs);
                end

				% Check if there is an RR interval larger than the 15% different than the prior interval 
                use_median = false;
        
                peakloc_diff = diff(peakloc);
                if length(peakloc_diff) < 50
                    pos_start = pos_start + interval_t;
                    disp('Warning: not enough peaks found')
                    continue;
                end
                first_peakloc=peakloc(1);
                if first_peakloc > 3000
                    disp('Warning: gap at the beginnning'); 
                    pos_start = pos_start + interval_t;
                    continue;
                end
                if use_median
                  median_PT=median(peakloc_diff);
                  upperlimit=median_PT*(1.0 + good_percent);
                  lowerlimit=median_PT*(1.0 - good_percent);
                  upperlimit_check=median_PT*(1.0 + check_percent);
                  lowerlimit_check=median_PT*(1.0 - check_percent);
                  upperlimit_check2=median_PT*(1.0 + check_percent2);
                  lowerlimit_check2=median_PT*(1.0 - check_percent2);
                  upperlimit_check3=median_PT*(1.0 + check_percent3);
                  lowerlimit_check3=median_PT*(1.0 - check_percent3);
                  upperlimit_check=median_PT*(1.0 + check_percent4);
                  lowerlimit_check4=median_PT*(1.0 - check_percent4);
                  upperlimit_check5=median_PT*(1.0 + check_percent5);
                  lowerlimit_check5=median_PT*(1.0 - check_percent5);
                  upperlimit_check6=median_PT*(1.0 + check_percent6);
                  lowerlimit_check6=median_PT*(1.0 - check_percent6);
                  upperlimit_check7=median_PT*(1.0 + check_percent7);
                  lowerlimit_check7=median_PT*(1.0 - check_percent7);
                  
                else
                  peakloc_diff_lag = vertcat([peakloc_diff(1)], peakloc_diff(1:length(peakloc_diff)-1));
                  upperlimit=peakloc_diff_lag*(1.0 + good_percent);
                  lowerlimit=peakloc_diff_lag*(1.0 - good_percent);
                  upperlimit_check=peakloc_diff_lag*(1.0 + check_percent);
                  lowerlimit_check=peakloc_diff_lag*(1.0 - check_percent);
                  upperlimit_check2=peakloc_diff_lag*(1.0 + check_percent2);
                  lowerlimit_check2=peakloc_diff_lag*(1.0 - check_percent2);
                  upperlimit_check3=peakloc_diff_lag*(1.0 + check_percent3);
                  lowerlimit_check3=peakloc_diff_lag*(1.0 - check_percent3);
                  upperlimit_check4=peakloc_diff_lag*(1.0 + check_percent4);
                  lowerlimit_check4=peakloc_diff_lag*(1.0 - check_percent4);
                  upperlimit_check5=peakloc_diff_lag*(1.0 + check_percent5);
                  lowerlimit_check5=peakloc_diff_lag*(1.0 - check_percent5);
                  upperlimit_check6=peakloc_diff_lag*(1.0 + check_percent6);
                  lowerlimit_check6=peakloc_diff_lag*(1.0 - check_percent6);
                  upperlimit_check7=peakloc_diff_lag*(1.0 + check_percent7);
                  lowerlimit_check7=peakloc_diff_lag*(1.0 - check_percent7);
                end

				[pos_d ~]=find(peakloc_diff>upperlimit | peakloc_diff<lowerlimit);
				[pos_d_check ~]=find(peakloc_diff>upperlimit_check | peakloc_diff<lowerlimit_check);
                [pos_d_check2 ~]=find(peakloc_diff>upperlimit_check2 | peakloc_diff<lowerlimit_check2);
				[pos_d_check3 ~]=find(peakloc_diff>upperlimit_check3 | peakloc_diff<lowerlimit_check3);
				[pos_d_check4 ~]=find(peakloc_diff>upperlimit_check | peakloc_diff<lowerlimit_check4);
                [pos_d_check5 ~]=find(peakloc_diff>upperlimit_check | peakloc_diff<lowerlimit_check5);
                [pos_d_check6 ~]=find(peakloc_diff>upperlimit_check | peakloc_diff<lowerlimit_check6);
                [pos_d_check7 ~]=find(peakloc_diff>upperlimit_check | peakloc_diff<lowerlimit_check7);
                
				% if good data is found, break the cycle (stop searching and exit the cycle)
				if isempty(pos_d)
                  break;
                else
                  disp('Num "bad" peaks:'), disp(length(pos_d))
                  % find largest index in pos_d
                  pos_d_last_index = max(pos_d);
                  % find corresponding index in peakloc and add to start position for next sliding window
                  if (length(peakloc) > pos_d_last_index)
                    % use the next peak as the new start location
                    next_valid_index_in_sequence = peakloc(pos_d_last_index+1);
                  else 
                    % skip to the end of the sequence
                    next_valid_index_in_sequence = interval_t;       
                  end
                  pos_start = pos_start + next_valid_index_in_sequence;
                  disp('Segment index of next window'), disp(next_valid_index_in_sequence)
                  disp('Length of data file:'), disp(sequence_length)
                  disp('Using algorithm:'), disp(algo_num)
                end
              end
              if isempty(pos_d) 
                % good data was found
                break; 
              else 
                % check if we have gone through each algorithm for rpeak detection 
                if algo_num < 3 %% change this to 3 when PCA method is working properly
                  algo_num = algo_num + 1;
                else 
                  break;
                end
              end
            end
			% Compute the values and write in the files if there is good data
			if isempty(pos_d)
              good_segment=start_index; 
              good_data=data;
              good_Rpeaks=peakloc;
              segments_to_check=pos_d_check;
              num_segments_to_check=length(pos_d_check);
              num_segments_to_check2=length(pos_d_check2);
              num_segments_to_check3=length(pos_d_check3);
              num_segments_to_check4=length(pos_d_check4);
              num_segments_to_check5=length(pos_d_check5);
              num_segments_to_check6=length(pos_d_check6);
              num_segments_to_check7=length(pos_d_check7);

             
              time_3min=((start_index-1.0)/interval_t)*3;
              Total_beats=length(peakloc);
				
				RRms_All=[];
				RRms=[];
				RRms=(peakloc/fs)*1000;
				RRms_All=diff(RRms);
	
			            % Basic parameters
			            meanRR=[]; sdRR=[]; RMSSD=[];
		        	    meanRR=mean(RRms_All(~isnan(RRms_All)));
			            sdRR=std(RRms_All(~isnan(RRms_All)));
			            RMSSD=sqrt((sum((diff(RRms_All)).^2))/(length(RRms_All)-1));


	                fid1 = fopen(excelName,'a');

	                fprintf(fid1,'%s \t',file_id);
	                fprintf(fid1,'%s \t',hour);
	                fprintf(fid1,'%d \t',time_3min);
	                fprintf(fid1,'%d \t',Total_beats);
	
	                fprintf(fid1,'%2.3f \t ',meanRR);
        	        fprintf(fid1,'%2.3f \t ',sdRR);
	                fprintf(fid1,'%2.3f \t ',RMSSD);
                    fprintf(fid1,'%d \t', algo_num);
	                fprintf(fid1,'%d \t ',num_segments_to_check);
                    fprintf(fid1,'%d \t ',num_segments_to_check2);
	                fprintf(fid1,'%d \t ',num_segments_to_check3);
	                fprintf(fid1,'%d \t ',num_segments_to_check4);
                    fprintf(fid1,'%d \t ',num_segments_to_check5);
	                fprintf(fid1,'%d \t ',num_segments_to_check6);
                    fprintf(fid1,'%d \t ',num_segments_to_check7);
                    
		            fprintf(fid1,'\n');
	        	    fclose(fid1);


			else

			        good_segment=[];
			        good_data=[];
			        good_Rpeaks=[];
                    segments_to_check=[];
                    num_segments_to_check=[];	
		                fid1 = fopen(excelName,'a');
	        	        fprintf(fid1,'%s \t',file_id);
	                	fprintf(fid1,'%s \t',hour);
	
		                fprintf(fid1,'\n');
		                fclose(fid1);

			end

			par_save(strcat(file_root,file_name),good_segment,good_data,good_Rpeaks,time_3min,Total_beats);
	
		end
	end
end


% function to save in the same matfile (this needs to be outside to be able to run in parallel
function par_save(filename,good_segment,good_data,good_Rpeaks,time_3min,Total_beats)
                save(filename,'-append','good_segment','good_data','good_Rpeaks','time_3min','Total_beats');
end

%Algorithms for detecting R peak (Pan tompkins, Principle component
%analysis, parabolic fitting) 
% R peak algo 1
function peakloc=rpeak_pan_tompkins(data,fs)
    %Rpeak detection using Pan-tomkins method 
    % threshold of 0.01 can be modified for each patient.
    peakloc=PeakDetection2(data,fs,[],[],[],0.01,[]);
    % Rcorrection based on R-peak neighborhood maximum to minimum
    th_sQRS=5;
    Rx_mxmn=[];

    for mnc=1:length(peakloc)
        Rx_mxmn(mnc,1)=max(data(max(1,peakloc(mnc)-th_sQRS):min(peakloc(mnc)+th_sQRS,length(data))));
        Rx_mxmn(mnc,2)=min(data(max(1,peakloc(mnc)-th_sQRS):min(peakloc(mnc)+th_sQRS,length(data))));
    end

    % R peak correction in case there is a missing peak
    Rx_loc=find((Rx_mxmn(:,1)-Rx_mxmn(:,2))>mean(Rx_mxmn(:,1)-Rx_mxmn(:,2))/2);
    Rx_temp=peakloc(Rx_loc);
    peakloc=Rx_temp;

    % R peak correction for absolute values
    peakloc = PeakCorrection(data,peakloc,'bestPeak');
end

% R peak algo 2
function Rx=rpeak_parabolic(data,fs)
    Rx=rpeak(data,fs);

    % Rcorrection based on R-peak neighborhood maximum to minimum
    th_sQRS=5;
    Rx_mxmn=[];
    for mnc=1:length(Rx)
        Rx_mxmn(mnc,1)=max(data(max(1,Rx(mnc)-th_sQRS):min(Rx(mnc)+th_sQRS,length(data))));
        Rx_mxmn(mnc,2)=min(data(max(1,Rx(mnc)-th_sQRS):min(Rx(mnc)+th_sQRS,length(data))));
    end

    Rx_loc=find((Rx_mxmn(:,1)-Rx_mxmn(:,2))>mean(Rx_mxmn(:,1)-Rx_mxmn(:,2))/2);
    Rx_temp=Rx(Rx_loc);
    Rx=Rx_temp;
    Rx = PeakCorrection(data,Rx,'bestPeak');
end

% R peak algo 3
function qrs_i_raw=rpeak_pca(data,fs)
    %offset = 100;
    
    [~,ECG_pca,~] = pca(data);
    ecg_avg = ECG_pca(:,1);

    if abs(prctile(ecg_avg(1:fix(length(ecg_avg)/2)),99))<abs(prctile(ecg_avg(1:fix(length(ecg_avg)/2)),1))
        ecg_avg = -ecg_avg;
    end

    % Take derivative to eccentuate R peaks before thresholding
    decg_avg = diff(ecg_avg);
    decg_avg2 = diff(decg_avg);
    thrsY = nanmean(decg_avg2)-1.5*std(decg_avg2);
    thrsX = nanmean(ecg_avg)+1.5*std(ecg_avg);

    %locate R peak
    ix = find(decg_avg(1:end-1)>0 & decg_avg(2:end)<0 & decg_avg2<thrsY & ecg_avg(3:end)>thrsX)+1;
    %ix = ix-offset;
    ix(ix<1) = [];

    qrs_i_raw=ix;
    e=3
end

