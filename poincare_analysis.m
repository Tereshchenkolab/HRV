function [SD1,SD2,SD1SD2ratio] = poincare_analysis(RR,num,flag)
%   This function computes standard deviations along the identity line and 
%   its perpendicular axis of the Poincare plot.
%   
%   RR: RR intervals in seconds.
%   num: the number of successive values for which the standard
%   deviations will be computed.  
%   For num = 0, the global standard deviations will be computed.
%   Otherwise SD1, SD2 and SD1SD2ratio are vectors of the same length as RR.
%   flag: 0 or 1 to specify normalization by n-1 or n. (default: 1)
%   

    if nargin<2 || isempty(num)
        num = 0;
    end    
    if nargin<3 || isempty(flag)
        flag = 1;
    end

    overlap = 1;
    
    RR = RR(:);

    X = [RR(1:end-1) RR(2:end)]';
    alpha = -45*pi/180;
    R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    XR = R*X;
    
    if num==0
        SD2 = nanstd(XR(1,:),flag,2);
        SD1 = nanstd(XR(2,:),flag,2);          
    else
        steps = ceil(num*(1-overlap));
        if steps>1
            j=1;
            ts1 = NaN(length(steps:steps:length(RR)-1),num);
            ts2 = ts1;            
            for i=steps:steps:length(RR)-1
                ts1(j,1:i-max(1,i-num+1)+1) = XR(1,max(1,i-num+1):i);
                ts2(j,1:i-max(1,i-num+1)+1) = XR(2,max(1,i-num+1):i);                
                j=j+1;
            end
            SD2_tmp = nanstd(ts1,flag,2);
            SD1_tmp = nanstd(ts2,flag,2); 
            
            SD2 = NaN(length(RR),1);
            SD1 = NaN(length(RR),1);
            SD2(steps+1:steps:length(RR)) = SD2_tmp;
            SD1(steps+1:steps:length(RR)) = SD1_tmp;
        else
            ts1 = NaN(length(RR),num);
            ts2 = NaN(length(RR),num);    
            for j=1:num
                ts1(j+1:end,j) = XR(1,1:end-j+1);
                ts2(j+1:end,j) = XR(2,1:end-j+1);        
            end
            SD2 = nanstd(ts1,flag,2);
            SD1 = nanstd(ts2,flag,2);
        end
    end
    SD1SD2ratio = SD1./SD2;
end