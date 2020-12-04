function [sf_1,sf_2,sf_3] = FiltroAlfa2(SIG,pos)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for j=1:72
        fs = 250;
        f_alfa=[7.99 8 13.98 13.99];                                           % Frequency Vector
        mags = [0 1 0];                                                 % Magnitude (Defines Passbands & Stopbands)
        devs = [0.0001 10^(0.1/20)-1 0.0001];                                            % Allowable Deviations
        [n,Wn,beta,ftype] = kaiserord(f_alfa,mags,devs,fs);              % Kaiser Window
        n = n + rem(n,2);
        hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'scale'); 
        p=(SIG(pos(j):pos(j)+1252,8)); 
        p2=(SIG(pos(j):pos(j)+1252,10)); 
        p3=(SIG(pos(j):pos(j)+1252,12));
        sf(:,j)=filter(hh,1,p);
        sf2(:,j)=filter(hh,1,p2);
        sf3(:,j)=filter(hh,1,p3);    
        sf_1=sf.';
        sf_2=sf2.';
        sf_3=sf3.'; 
    end
end