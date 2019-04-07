function hrv_rmssd = RMSSDcalc(RR,num,flag)  
%RMSSDcalc Root Mean Square of Successive Differences.
%   RR: RR intervals in seconds.
%   num: the number of successive values for which the standard
%   deviations will be computed.  
%   For num = 0, the global standard deviations will be computed.
%   Otherwise hrv_rmssd is vector of the same length as RR.
%   flag: 0 or 1 to specify normalization by n-1 or n. (default: 1)
%   

    if nargin<2 || isempty(num)
        num = 0;
    end
    if nargin<3 || isempty(flag)
        flag = 1; %The flag is 0 or 1 to specify normalization by n-1 or n.
    end 
    overlap = 1;
    
    RR = RR(:);

    dRR = diff(RR).^2;
    if num==0
        hrv_rmssd = sqrt(nansum(dRR)./(sum(~isnan(dRR))-1+flag));
    else
        if ceil(num*(1-overlap))>1
            j=1;
            ts = NaN(length(ceil(num*(1-overlap)):ceil(num*(1-overlap)):length(dRR)),num);
            for i=ceil(num*(1-overlap)):ceil(num*(1-overlap)):length(dRR)
                ts(j,1:(1+i-max(1,(i-num+1)))) = dRR(max(1,(i-num+1)):i);
                j=j+1;
            end 
            samplesize = sum(~isnan(ts),2);
            hrv_rmssd_tmp = sqrt(nansum(ts,2)./(samplesize-1+flag)); 
            hrv_rmssd_tmp(samplesize<5) = NaN;
            
            hrv_rmssd = NaN(length(RR),1);  
            hrv_rmssd(ceil(num*(1-overlap))+1:ceil(num*(1-overlap)):length(RR)) = hrv_rmssd_tmp;  
        else
            ts = NaN(length(RR),num);
            for j=1:num
                ts(j+1:end,j) = dRR(1:end-j+1);
            end    
            samplesize = sum(~isnan(ts),2);
            hrv_rmssd = sqrt(nansum(ts,2)./(samplesize-1+flag));
            hrv_rmssd(samplesize<5) = NaN;
        end
    end
end