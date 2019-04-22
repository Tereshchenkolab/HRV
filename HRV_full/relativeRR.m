function [r,rrHRV]=relativeRR(RRINT)

r=[];
rrHRV=[];
rr=[];

% Relative RR-intervals (Marcus Vollmer's method)
j1=1;
for j=1:length(RRINT)-1
    r(j)= 2*(RRINT(j+1)-RRINT(j))/(RRINT(j+1)+RRINT(j));
    if abs(r(j))<0.2
        rr(j1)=r(j);
        j1 = j1+1;
    end
end
%     figure
%     plot(rr(1:end-1),rr(2:end),'o'); 
%     title('Return map of rr');

% Center point
if ~isempty(rr)
    c = [mean(rr(1:end-1)),mean(rr(2:end))];

    dd=[];
    for j=1:length(rr)-1
        dd(j) = norm([rr(j+1),rr(j)]-c);
    end
    
    if ~isempty(dd)
        rrHRV = median(dd(~isnan(dd)));
    end
end

