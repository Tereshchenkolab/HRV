function[SDmk_Two,allSym_Two]=calc_SymDynMK_WLtwo(RRms)

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

% coding into 2 letter word type probabilities
R0=0;
RI=0;
RD=0;
for i=1:length(symRR)-1
    word=symRR(i:i+1);
    if (word(1)==word(2)), R0=R0+1; end
    if (word(1)>word(2)), RI=RI+1; end
    if (word(1)<word(2)), RD=RD+1; end
end

SDmk_Two=[];
allSym_Two=[];
allSym_Two=R0+RI+RD;
% SDmk_Two.R0=R0/allSym_Two*100;
% SDmk_Two.RI=RI/allSym_Two*100;
% SDmk_Two.RD=RD/allSym_Two*100;
SDmk_Two.R0=R0;
SDmk_Two.RI=RI;
SDmk_Two.RD=RD;
