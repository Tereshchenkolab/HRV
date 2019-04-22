function[SDmk,RRalt_mag,allSym]=calc_SymDynMK(RRms)

symRR=[];
for i=1:length(RRms)-1
    if RRms(i+1)-RRms(i)<0
        symRR(i)=-1;
    elseif RRms(i+1)-RRms(i)>0
        symRR(i)=1;
    else
        symRR(i)=0;
    end
end

%symRR as string
symRRstr='';
for i=1:length(symRR)
    symRRstr=strcat(symRRstr,num2str(symRR(i)));
end

% coding into 3 letter word type probabilities
R0=0;
RL=0;
RUL=0;
R1=0;
ALT=0;
RRalt_mag=[];
for i=1:length(symRR)-2
    word=symRR(i:i+2);
    if (word(1)==word(2))&& (word(2)==word(3)), R0=R0+1; end
    if ((word(1)>word(2))&&(word(2)>word(3))) || ((word(1)<word(2))&& (word(2)<word(3))), RL=RL+1; end
    if ((word(1)<word(2))&&(word(2)>word(3))) || ((word(1)>word(2)) &&(word(2)<word(3))), RUL=RUL+1; end
    if ((word(1)==word(2))&&(word(2)~=word(3))) || ((word(1)~=word(2))&&(word(2)==word(3))), R1=R1+1; end
    if max(ismember(word,0))==0
        if ((word(1)<word(2))&&(word(2)>word(3))) || ((word(1)>word(2)) &&(word(2)<word(3))) 
            ALT=ALT+1;
            RRalt_mag(ALT)=(abs(RRms(i+2)-RRms(i+1))+abs(RRms(i+3)-RRms(i+2)))/2;
        end
    end
end

SDmk=[];
allSym=[];
allSym=R0+R1+RUL+RL;
% SDmk.R0=R0/allSym*100;
% SDmk.R1=R1/allSym*100;
% SDmk.LR2=RL/allSym*100;
% SDmk.ULR2=RUL/allSym*100;
% SDmk.ALT=ALT/allSym*100;
SDmk.R0=R0;
SDmk.R1=R1;
SDmk.LR2=RL;
SDmk.ULR2=RUL;
SDmk.ALT=ALT;