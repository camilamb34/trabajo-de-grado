function [Delta] = FiltroDelta(SIG,a)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for j=1:72
        fs = 250;
        f_delta=[7 8 40 41];                                           % Frequency Vector
        mags = [0 1 0];                                                 % Magnitude (Defines Passbands & Stopbands)
        devs = [0.0001 10^(0.1/20)-1 0.0001];                                            % Allowable Deviations
        [n,Wn,beta,ftype] = kaiserord(f_delta,mags,devs,fs);              % Kaiser Window
        n = n + rem(n,2);
        hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'scale'); 
        p=(SIG(a(j):a(j)+1252,8)); 
        p2=(SIG(a(j):a(j)+1252,10)); 
        p3=(SIG(a(j):a(j)+1252,12));
        sf=filter(hh,1,p);
        sf2=filter(hh,1,p2);
        sf3=filter(hh,1,p3);    
        Delta(:,j) = (sf+sf2+sf3)/3;
        TransDelta=(Delta).';
    end
end

