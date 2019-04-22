function[SDPorta,all]=calc_SymDynPorta(rr)

lowbound=min(rr);
upbound=max(rr);

stepwidth=(upbound-lowbound)/6;
% coding into 6 different symbols
symrr=[];
for i=1:length(rr)
    if rr(i) >= lowbound && rr(i)< lowbound+stepwidth, symrr(i)=0;end
    if rr(i) >= lowbound+stepwidth && rr(i)< lowbound+stepwidth*2, symrr(i)=1;end
    if rr(i) >= lowbound+stepwidth*2 && rr(i)< lowbound+stepwidth*3, symrr(i)=2;end 
    if rr(i) >= lowbound+stepwidth*3 && rr(i)< lowbound+stepwidth*4, symrr(i)=3;end 
    if rr(i) >= lowbound+stepwidth*4 && rr(i)< lowbound+stepwidth*5, symrr(i)=4;end 
    if rr(i) >= lowbound+stepwidth*5 && rr(i)<= upbound, symrr(i)=5;end
end

% coding into 3 letter word type probabilities
V0=0;
VL=0;
VUL=0;
V1=0;
for i=1: length(rr)-2
    word=symrr(i:i+2);
    if (word(1)==word(2))&& (word(2)==word(3)), V0=V0+1; end
    if ((word(1)>word(2))&&(word(2)>word(3))) || ((word(1)<word(2))&& (word(2)<word(3))), VL=VL+1; end
    if ((word(1)<word(2))&&(word(2)>word(3))) || ((word(1)>word(2)) &&(word(2)<word(3))), VUL=VUL+1; end
    if ((word(1)==word(2))&&(word(2)~=word(3))) || ((word(1)~=word(2))&&(word(2)==word(3))), V1=V1+1; end
end

SDPorta=[];
all=[];
all=V0+V1+VUL+VL;
% SDPorta.V0=V0/all*100;
% SDPorta.V1=V1/all*100;
% SDPorta.LV2=VL/all*100;
% SDPorta.ULV2=VUL/all*100;
% SDPorta.ALT=ALT/all*100;
SDPorta.V0=V0;
SDPorta.V1=V1;
SDPorta.LV2=VL;
SDPorta.ULV2=VUL;


