function Peak_corr = PeakCorrection(signal,Peak,method) 
    vicinity1 = 5; %14
%     vicinity2 = 135;
    Peak_corr = [];
%     QONx      = [mFile.q_points_x,mFile.q_points_y,mFile.q_points_z];
    r         = [];
    Qcorrx    = [];
    Qcorry    = [];
    Qcorrz    = [];
if strcmp(method, 'bestPeak')
  % This method will pick either all positive peaks or all negative peaks depending on
  %  which results in the largest average manhattan (L1) distance from the median 
  posPeaks=PeakCorrection(signal,Peak,'posPeak');
  negPeaks=PeakCorrection(signal,Peak,'negPeak');
  
  if negPeaks(1) < 1
    negPeaks=negPeaks(2:length(negPeaks));    
  end
  if posPeaks(1) < 1
    posPeaks=posPeaks(2:length(posPeaks));    
  end
  
  % compute the median across the entire signal
  signalMedian=median(signal);
  % get the signal values corresponding to the indices of the positive and negative peaks
  posPeakValues=signal(posPeaks);
  negPeakValues=signal(negPeaks);
  % compute L1 distance from the mean
  posDeviation=mean(abs(posPeakValues-signalMedian));
  negDeviation=mean(abs(negPeakValues-signalMedian));
  % return whichever resulted in the largest total deviation
  if posDeviation >= negDeviation
    Peak_corr=posPeaks;
  else
    Peak_corr=negPeaks;    
  end
  
elseif strcmp(method,'absPeak')    
    if size(signal,2)==1
        for i=1:length(Peak)
            if Peak(i)<=vicinity1
                Peak_corr=Peak(i);
                
            else
%                 if i==length(Peak)
%                     window1       = Peak(i)-vicinity1:length(signal);
%                 else
                    window1       = Peak(i)-vicinity1:Peak(i)+vicinity1;
%                 end
    %             window1   = Peak(i)-vicinity1:Peak(i);
    %             [~,mm1] = max(signal(window1,1));

                window1(window1<=0)=[];
                window1(window1>=length(signal))=[];                
                [~,Peak_temp] = max(abs(signal(window1)));
                Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];
            end

    % %% QRS onset
    %         if window1(mm1)~=Peak(i,1)
    %             [~,Qtemp] = min(signal(window1(mm1)-vicinity2:window1(mm1),1));
    %             Qcorrx     = [Qcorrx;window1(mm1)-vicinity2+Qtemp-1]; 
    %         else
    %             [~,Qtemp] = min(signal(window1,1));
    %             Qcorrx     = [Qcorrx;Peak(i,1)-vicinity2+Qtemp-1]; 
    %         end 
        end
    else
    
        for i=1:length(Peak)            
            window1   = Peak(i)-vicinity1:Peak(i)+vicinity1;
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
            [~,Peak_temp] = max(abs(signal(window1)));
            Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];            
            window2   = Peak(i,2)-vicinity2:Peak(i,2);
            window2(window2<=0)=[];
            window2(window2>=length(signal))=[];
            [~,mm2] = max(signal(window2,2));
            if window2(mm2)~=Peak(i,2)
                [~,Qtemp] = min(signal(window2(mm2)-vicinity2:window2(mm2),2));
                Qcorry     = [Qcorry;window2(mm2)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window2,2));
                Qcorry     = [Qcorry;Peak(i,2)-vicinity2+Qtemp-1]; 
            end 

            window3   = Peak(i,3)-vicinity2:Peak(i,3);
            window3(window3<=0)=[];
            window3(window3>=length(signal))=[];
            [~,mm3] = max(signal(window3,3));            
            if window3(mm3)~=Peak(i,3)
                [~,Qtemp] = min(signal(window3(mm3)-vicinity2:window3(mm3),3));
                Qcorrz     = [Qcorrz;window3(mm3)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window3,3));
                Qcorrz     = [Qcorrz;Peak(i,3)-vicinity2+Qtemp-1]; 
            end 
        end
        q_points_x = Qcorrx;
        q_points_y = Qcorry;
        q_points_z = Qcorrz;
    end
elseif strcmp(method,'negPeak')
    if size(signal,2)==1
        for i=1:length(Peak)
            if i==length(Peak) && Peak(i)+vicinity1>length(Peak)
                window1       = Peak(i)-vicinity1:length(signal);
            else
                window1       = Peak(i)-vicinity1:Peak(i)+vicinity1;
            end
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
            [~,Peak_temp] = min(signal(window1));
            Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];
    % %% QRS onset
    %         if window1(mm1)~=Peak(i,1)
    %             [~,Qtemp] = min(signal(window1(mm1)-vicinity2:window1(mm1),1));
    %             Qcorrx     = [Qcorrx;window1(mm1)-vicinity2+Qtemp-1]; 
    %         else
    %             [~,Qtemp] = min(signal(window1,1));
    %             Qcorrx     = [Qcorrx;Peak(i,1)-vicinity2+Qtemp-1]; 
    %         end 
        end
    else
    
        for i=1:length(Peak)            
            window1   = Peak(i)-vicinity1:Peak(i)+vicinity1;
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
            [~,Peak_temp] = min(signal(window1));
            Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];            
            window2   = Peak(i,2)-vicinity2:Peak(i,2);
            window2(window2<=0)=[];
            window2(window2>=length(signal))=[];
            [~,mm2] = max(signal(window2,2));
            if window2(mm2)~=Peak(i,2)
                [~,Qtemp] = min(signal(window2(mm2)-vicinity2:window2(mm2),2));
                Qcorry     = [Qcorry;window2(mm2)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window2,2));
                Qcorry     = [Qcorry;Peak(i,2)-vicinity2+Qtemp-1]; 
            end 

            window3   = Peak(i,3)-vicinity2:Peak(i,3);
            window3(window3<=0)=[];
            window3(window3>=length(signal))=[];
            [~,mm3] = max(signal(window3,3));            
            if window3(mm3)~=Peak(i,3)
                [~,Qtemp] = min(signal(window3(mm3)-vicinity2:window3(mm3),3));
                Qcorrz     = [Qcorrz;window3(mm3)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window3,3));
                Qcorrz     = [Qcorrz;Peak(i,3)-vicinity2+Qtemp-1]; 
            end 
        end
        q_points_x = Qcorrx;
        q_points_y = Qcorry;
        q_points_z = Qcorrz;
    end    
elseif strcmp(method,'posPeak')
    if size(signal,2)==1
        for i=1:length(Peak)
            window1   = Peak(i)-vicinity1:Peak(i);
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
window1
            [~,mm1] = max(signal(window1,1));

            window1       = Peak(i)-vicinity1:Peak(i)+vicinity1;
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
            [~,Peak_temp] = max(signal(window1));
            Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];
    % %% QRS onset
    %         if window1(mm1)~=Peak(i,1)
    %             [~,Qtemp] = min(signal(window1(mm1)-vicinity2:window1(mm1),1));
    %             Qcorrx     = [Qcorrx;window1(mm1)-vicinity2+Qtemp-1]; 
    %         else
    %             [~,Qtemp] = min(signal(window1,1));
    %             Qcorrx     = [Qcorrx;Peak(i,1)-vicinity2+Qtemp-1]; 
    %         end 
        end
    else
    
        for i=1:length(Peak)            
            window1   = Peak(i)-vicinity1:Peak(i)+vicinity1;
            window1(window1<=0)=[];
            window1(window1>=length(signal))=[];
            [~,Peak_temp] = max(signal(window1));
            Peak_corr     = [Peak_corr;Peak(i)-vicinity1+Peak_temp-1];            
            window2   = Peak(i,2)-vicinity2:Peak(i,2);
            window2(window2<=0)=[];
            window2(window2>=length(signal))=[];
            [~,mm2] = max(signal(window2,2));
            if window2(mm2)~=Peak(i,2)
                [~,Qtemp] = min(signal(window2(mm2)-vicinity2:window2(mm2),2));
                Qcorry     = [Qcorry;window2(mm2)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window2,2));
                Qcorry     = [Qcorry;Peak(i,2)-vicinity2+Qtemp-1]; 
            end 

            window3   = Peak(i,3)-vicinity2:Peak(i,3);
            window3(window3<=0)=[];
            window3(window3>=length(signal))=[];
            [~,mm3] = max(signal(window3,3));            
            if window3(mm3)~=Peak(i,3)
                [~,Qtemp] = min(signal(window3(mm3)-vicinity2:window3(mm3),3));
                Qcorrz     = [Qcorrz;window3(mm3)-vicinity2+Qtemp-1]; 
            else
                [~,Qtemp] = min(signal(window3,3));
                Qcorrz     = [Qcorrz;Peak(i,3)-vicinity2+Qtemp-1]; 
            end 
        end
        q_points_x = Qcorrx;
        q_points_y = Qcorry;
        q_points_z = Qcorrz;
    end    
    
end




    
end
