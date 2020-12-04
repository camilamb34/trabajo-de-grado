function [sf_1,sf_2,sf_3] = FiltroBeta2(SIG,pos)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for j=1:72
        fs = 250;
        f_beta=[13.99 14 30 30.01];                                           % Frequency Vector
        mags = [0 1 0];                                                 % Magnitude (Defines Passbands & Stopbands)
        devs = [0.0001 10^(0.1/20)-1 0.0001];                                            % Allowable Deviations
        [n,Wn,beta,ftype] = kaiserord(f_beta,mags,devs,fs);              % Kaiser Window
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