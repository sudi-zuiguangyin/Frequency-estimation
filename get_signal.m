function [sig, ref_pha] = get_signal(para, t, dist)
fundamp = para(1);
amp_mod = para(2);
fre_mod = para(3);
iniang_mod = para(4);
deltaf = para(5);
ang_mod = para(6);
dfre = para(7);
iniang = para(8);

%get sampling values
sig = sqrt(2)*fundamp*(1+amp_mod*cos(2*pi*fre_mod*t+iniang_mod)) .* cos(2*pi*50*t + 2*pi*deltaf*t + ...
    ang_mod*cos(2*pi*fre_mod*t-pi+iniang_mod) + pi*dfre*t.*t + iniang);%test signal model
%calculate reference values
amp = fundamp*(1+amp_mod*cos(2*pi*fre_mod*t+iniang_mod));
ang = 2*pi*deltaf*t + ang_mod*cos(2*pi*fre_mod*t-pi+iniang_mod) + pi*dfre*t.*t + iniang;
ang = ang*180/pi;
ang = ang - round(ang/360)*360;
fre = 50 + deltaf - fre_mod*ang_mod*sin(2*pi*fre_mod*t-pi+iniang_mod) + dfre*t;
corof = -2*pi*fre_mod*fre_mod*ang_mod*cos(2*pi*fre_mod*t-pi+iniang_mod) + dfre;

amp = amp(1: dist: length(amp));
ang = ang(1: dist: length(ang));
fre = fre(1: dist: length(fre));
corof = corof(1: dist: length(corof));
ref_pha = [amp', ang', fre', corof'];












