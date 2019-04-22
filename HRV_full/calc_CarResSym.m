function[SDmk,allSym]=calc_CarResSym(RRms,ph_X_R)

% Symbolics for RR
% +1 means RR increasing, -1 means RR decreasing, 0 means RR constant
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


% Symbolics for Respiratory phases
% +1 means inspiration, -1 means expiration, 0 means transition
symRsPh=[]; 
for i=1:length(ph_X_R)-1
    if ph_X_R(i)<0 && ph_X_R(i+1)<0
        symRsPh(i)=1;
    elseif ph_X_R(i)>0 && ph_X_R(i+1)>0
        symRsPh(i)=-1;
    else
        symRsPh(i)=0;
    end
end


% coding into 3 letter word type probabilities
R0=0;
RL=0;
RUL=0;
R1=0;
ALT=0;
ALT_exResTr=0; %excluding respiratory phase transition

for i=1:length(symRR)-2
    word=symRR(i:i+2);
    wordRsPh=symRsPh(i:i+2);
    if (word(1)==word(2))&& (word(2)==word(3)), R0=R0+1; end
    if ((word(1)>word(2))&&(word(2)>word(3))) || ((word(1)<word(2))&& (word(2)<word(3))), RL=RL+1; end
    if ((word(1)<word(2))&&(word(2)>word(3))) || ((word(1)>word(2)) &&(word(2)<word(3))), RUL=RUL+1; end
    if ((word(1)==word(2))&&(word(2)~=word(3))) || ((word(1)~=word(2))&&(word(2)==word(3))), R1=R1+1; end
    if max(ismember(word,0))==0
        if ((word(1)<word(2))&&(word(2)>word(3))) || ((word(1)>word(2)) &&(word(2)<word(3))) 
            ALT=ALT+1;
            if max(ismember(wordRsPh,0))==0
                ALT_exResTr=ALT_exResTr+1;
            end
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
SDmk.ALT_exResTr=ALT_exResTr;