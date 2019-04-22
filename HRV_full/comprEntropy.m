function Hc = comprEntropy(RRms,wl,lb)

%wl: sliding window length
%lb: lookahead buffer length
Hc=zeros(wl,lb);
pointer=1;

for w=1:1:wl    
    for b=1:1:lb
        Matr=[]; matCT=1; %matrix
        if w+b<=length(RRms)
            pointer=w+1;
            rrCT=pointer-1+b;
            brCT=1; %Break count
            while rrCT<=length(RRms) && brCT<=30
                VN=[];
                for v=1:w
                    VN(v)=NaN;
                    for n=1:b
                        if RRms(pointer+n-1)==RRms(pointer-w+(v-1)+(n-1))
                            continue;
                        else
                            VN(v)=n-1;
                            break;
                        end
                    end
                    if isnan(VN(v))
                        VN(v)=n;
                    end
                end
                [mxN, mxV]=max(VN);
                Matr(matCT,1)=mxV-1;
                Matr(matCT,2)=mxN;
                if pointer+mxN<=length(RRms)
                    Matr(matCT,3)=RRms(pointer+mxN);
                else
                    Matr(matCT,3)=0;
                end
                matCT=matCT+1;
                
                pointer=pointer+mxN+1;
                rrCT=pointer-1+b;
                
                brCT=brCT+1;
            end
        end
        if ~isempty(Matr)
            Hc(w,b)=length(Matr(:,1))/length(RRms);
        end
    end
end