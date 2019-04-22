function[HRV_SymDyn]=calc_HRV_SymDyn_Javorka(tacho)

        samples = length(tacho);
        diff=0;
        amean(1:samples) = mean(tacho); 
        a = 0.05; %0.1*std(tacho);%
        SD=std(tacho);
        rrdiff=tacho(2:end)-tacho(1:end-1);
        rrdiff=rrdiff.^2;
        rmssd=sqrt(mean(rrdiff));
        %%%%%% Standard %%%%%%%%%%%
        % settings_nld(17) = 0.05
        % settings_nld(18) = 1.5
        % settings_nld(19) = 0.01
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        limit1 = rmssd;     % plvar2 phvar2
        limit2 = rmssd;     % plvar5 phvar5
        limit3 = 10;    % plvar10 phvar10
        limit4 = 20;    % plvar20 phvar20
        
        %%%%%%%%%%%%%% variance of differences
        diff = tacho(2:samples)-tacho(1:samples-1);
        samples = length(diff);
        diffvar = std(diff);
    
       counter=1;
   
       symbols = zeros(1,samples-1);  % Symbolvektor1
       symbols2 = zeros(1,samples-1); % Symbolvektor2 - Suche nach "000000" bzw "111111" für plvar phvar
       symbols3 = zeros(1,samples-1);
       symbols4 = zeros(1,samples-1);
       symbols5 = zeros(1,samples-1);

          for i=1:samples-1
              if amean(i) < tacho(i) && tacho(i) <= SD+amean(i), symbols(i)=0; end
              if SD+amean(i) < tacho(i), symbols(i)=1; end
              if amean(i)-SD < tacho(i) && tacho(i) <= amean(i), symbols(i)=2; end
              if tacho(i) <= amean(i)-SD, symbols(i)=3; end
              if abs(tacho(i+1)-tacho(i)) < limit1, symbols2(i)=0; end
              if abs(tacho(i+1)-tacho(i)) >= limit1, symbols2(i)=1; end
              if abs(tacho(i+1)-tacho(i)) < limit2, symbols3(i)=0; end
              if abs(tacho(i+1)-tacho(i)) >= limit2, symbols3(i)=1; end
              if abs(tacho(i+1)-tacho(i)) < limit3, symbols4(i)=0; end
              if abs(tacho(i+1)-tacho(i)) >= limit3, symbols4(i)=1; end
              if abs(tacho(i+1)-tacho(i)) < limit4, symbols5(i)=0; end
              if abs(tacho(i+1)-tacho(i)) >= limit4, symbols5(i)=1; end
          end

     
%     words = zeros(1,64);
%     muster = ['000'; '001'; '002'; '003'; '010'; '011'; '012'; '013'; '020'; '021'; '022'; '023'; '030'; '031'; '032'; '033'; '100'; '101'; '102'; '103'; '110'; '111'; '112'; '113'; '120'; '121'; '122'; '123'; '130'; '131'; '132'; '133'; '200'; '201'; '202'; '203'; '210'; '211'; '212'; '213'; '220'; '221'; '222'; '223'; '230'; '231'; '232'; '233'; '300'; '301'; '302'; '303'; '310'; '311'; '312'; '313'; '320'; '321'; '322'; '323'; '330'; '331'; '332'; '333'];
%     symbolstring = zeros([1 length(symbols)]);
%     symbolstring = num2str(symbols(:));
%     wordsymbol = zeros([1 length(symbolstring)-2]);
%     wordcounter=1;
%     wpsum02(counter)=0;
%     wpsum13(counter)=0;
%     wsdvar(counter)=0; % 02/2002
%     ShannEntr(counter)=0; % 02/2002
%     Renyi025(counter)=0;
%     Renyi4(counter)=0;
%     for r=1:length(symbolstring)-2
%         symbmuster=strcat(symbolstring(r),symbolstring(r+1),symbolstring(r+2));
%         if length(find(symbmuster=='0' | symbmuster=='2' ))==3
%            wpsum02(counter)=wpsum02(counter)+1; 
%        end
%        if length(find(symbmuster=='1' | symbmuster=='3' ))==3
%           wpsum13(counter)=wpsum13(counter)+1; 
%        end
%        for i=1:64
%            muster2 = strcat(muster(i,1),muster(i,2),muster(i,3));     
%            if symbmuster==muster2
%               words(i)=words(i)+1; 
%            end
%        end
%     end
    words = zeros(1,64);
    muster = [0 0 0; 0 0 1; 0 0 2; 0 0 3; 0 1 0; 0 1 1; 0 1 2; 0 1 3; 0 2 0; 0 2 1; 0 2 2; 0 2 3; 0 3 0; 0 3 1; 0 3 2; 0 3 3; 1 0 0; 1 0 1; 1 0 2; 1 0 3; 1 1 0; 1 1 1; 1 1 2; 1 1 3; 1 2 0; 1 2 1; 1 2 2; 1 2 3; 1 3 0; 1 3 1; 1 3 2; 1 3 3; 2 0 0; 2 0 1; 2 0 2; 2 0 3; 2 1 0; 2 1 1; 2 1 2; 2 1 3; 2 2 0; 2 2 1; 2 2 2; 2 2 3; 2 3 0; 2 3 1; 2 3 2; 2 3 3; 3 0 0; 3 0 1; 3 0 2; 3 0 3; 3 1 0; 3 1 1; 3 1 2; 3 1 3; 3 2 0; 3 2 1; 3 2 2; 3 2 3; 3 3 0; 3 3 1; 3 3 2; 3 3 3];
    wordsymbol = zeros([1 length(symbols)-2]);
    wordcounter=1;
    wpsum02(counter)=0;
    wpsum13(counter)=0;
    wsdvar(counter)=0; % 02/2002
    ShannEntr(counter)=0; % 02/2002
    Renyi025(counter)=0;
    Renyi4(counter)=0;
    plvar(counter)=0;
    phvar(counter)=0;
    
%     for r=1:length(symbolstring)-6   
%         symbmuster=strcat(symbolstring(r),symbolstring(r+1),symbolstring(r+2));
%         eins=find(symbmuster=='1');
%         drei=find(symbmuster=='3');
%         if not(isempty(eins)) & not(isempty(drei)) 
%             if eins(1)<drei(1)      % im Wort der Länge 3 kommt eine 1 vor einer drei
%                 leins = length(eins) + length(drei); % 02/2002
%                 if leins==3, wordsymbol(wordcounter)=3; wordcounter=wordcounter+1; end
%                 if leins==2, wordsymbol(wordcounter)=2; wordcounter=wordcounter+1; end
%                 if leins==1, wordsymbol(wordcounter)=1; wordcounter=wordcounter+1; end  
%             end
%             if eins(1)>drei(1) 
%                 ldrei = length(drei) + length(eins); % 02/2002
%                 if ldrei==3, wordsymbol(wordcounter)=-3; wordcounter=wordcounter+1; end
%                 if ldrei==2, wordsymbol(wordcounter)=-2; wordcounter=wordcounter+1; end
%                 if ldrei==1, wordsymbol(wordcounter)=-1; wordcounter=wordcounter+1; end
%             end
%         elseif not(isempty(eins)) && isempty(drei)
%                 leins = length(eins);
%                 if leins==3, wordsymbol(wordcounter)=3; wordcounter=wordcounter+1; end
%                 if leins==2, wordsymbol(wordcounter)=2; wordcounter=wordcounter+1; end
%                 if leins==1, wordsymbol(wordcounter)=1; wordcounter=wordcounter+1; end
%         elseif isempty(eins) && not(isempty(drei))
%                 ldrei = length(drei);
%                 if ldrei==3, wordsymbol(wordcounter)=-3; wordcounter=wordcounter+1; end
%                 if ldrei==2, wordsymbol(wordcounter)=-2; wordcounter=wordcounter+1; end
%                 if ldrei==1, wordsymbol(wordcounter)=-1; wordcounter=wordcounter+1; end
%         else
%                 wordsymbol(wordcounter)=0; wordcounter=wordcounter+1;
%         end 
%     end
%********* neu optimiert *****************************************************************
    for r=1:length(symbols)-2
        if length(find(symbols(r:r+2)==0 | symbols(r:r+2)==2 ))==3
           wpsum02(counter)=wpsum02(counter)+1; 
       end
       if length(find(symbols(r:r+2)==1 | symbols(r:r+2)==3 ))==3
          wpsum13(counter)=wpsum13(counter)+1; 
       end
       for i=1:64     
           if symbols(r:r+2)==muster(i,:)
              words(i)=words(i)+1; break;
           end
       end
        eins=find(symbols(r:r+2)==1);
        drei=find(symbols(r:r+2)==3);
        if not(isempty(eins)) & not(isempty(drei)) 
            if eins(1)<drei(1)      % im Wort der Länge 3 kommt eine 1 vor einer drei
                leins = length(eins) + length(drei); % 02/2002
                if leins==3, wordsymbol(wordcounter)=3; wordcounter=wordcounter+1; end
                if leins==2, wordsymbol(wordcounter)=2; wordcounter=wordcounter+1; end
                if leins==1, wordsymbol(wordcounter)=1; wordcounter=wordcounter+1; end  
            end
            if eins(1)>drei(1) 
                ldrei = length(drei) + length(eins); % 02/2002
                if ldrei==3, wordsymbol(wordcounter)=-3; wordcounter=wordcounter+1; end
                if ldrei==2, wordsymbol(wordcounter)=-2; wordcounter=wordcounter+1; end
                if ldrei==1, wordsymbol(wordcounter)=-1; wordcounter=wordcounter+1; end
            end
        elseif not(isempty(eins)) & isempty(drei)
                leins = length(eins);
                if leins==3, wordsymbol(wordcounter)=3; wordcounter=wordcounter+1; end
                if leins==2, wordsymbol(wordcounter)=2; wordcounter=wordcounter+1; end
                if leins==1, wordsymbol(wordcounter)=1; wordcounter=wordcounter+1; end
        elseif isempty(eins) & not(isempty(drei))
                ldrei = length(drei);
                if ldrei==3, wordsymbol(wordcounter)=-3; wordcounter=wordcounter+1; end
                if ldrei==2, wordsymbol(wordcounter)=-2; wordcounter=wordcounter+1; end
                if ldrei==1, wordsymbol(wordcounter)=-1; wordcounter=wordcounter+1; end
        else
                wordsymbol(wordcounter)=0; wordcounter=wordcounter+1;
        end 
    end
    words = words./sum(words);
    ForbWord(counter) = 0;
    for i=1:64
        if words(i) == 0
           ShannEntr(counter)=ShannEntr(counter); % 02/2002
           Renyi025(counter)=Renyi025(counter);
           Renyi4(counter)=Renyi4(counter);
        else
           ShannEntr(counter)=(words(i)*log(words(i))) + ShannEntr(counter); % 02/2002
           Renyi025(counter)=words(i)^0.25 + Renyi025(counter);
           Renyi4(counter)=words(i)^4 + Renyi4(counter);
        end
        if words(i) <= 0.001
           ForbWord(counter) = ForbWord(counter)+1; 
        end
    end
    ShannEntr(counter)= -1 * ShannEntr(counter);
    Renyi025(counter)=4/3*log(Renyi025(counter));
    Renyi4(counter)=-1/3*log(Renyi4(counter));
    wpsum02(counter) = wpsum02(counter)/samples;
    wpsum13(counter) = wpsum13(counter)/samples;
    wsdvar(counter) = std(wordsymbol);
    counter=counter+1;
%     plvar2=find6(symbols2,0)/samples;    
%     phvar2=find6(symbols2,1)/samples;  
    plvar5=find6(symbols3,0)/samples;    
    phvar5=find6(symbols3,1)/samples;
%     plvar10=find6(symbols4,0)/samples;    
%     phvar10=find6(symbols4,1)/samples;
%     plvar20=find6(symbols5,0)/samples;    
%     phvar20=find6(symbols5,1)/samples;

    %--------------------------------------------------------------------------------------------------------------------------------
    % create output variable
    HRV_SymDyn.ShannonEntr_SD=ShannEntr; 
    HRV_SymDyn.ForbWord_SD=ForbWord;
    HRV_SymDyn.wpsum02_SD=wpsum02;
    HRV_SymDyn.wpsum13_SD=wpsum13;          
    HRV_SymDyn.wsdvar_SD=wsdvar;          
%     HRV_SymDyn.plvar2=plvar2;           
%     HRV_SymDyn.phvar2=phvar2;           
    HRV_SymDyn.plvarRMSSD=plvar5;           
    HRV_SymDyn.phvarRMSSD=phvar5;          
%     HRV_SymDyn.plvar10=plvar10;         
%     HRV_SymDyn.phvar10=phvar10;       
%     HRV_SymDyn.plvar20=plvar20;       
%     HRV_SymDyn.phvar20=phvar20;        
    HRV_SymDyn.Renyi025_SD=Renyi025;       
    HRV_SymDyn.Renyi4_SD=Renyi4;         
 


