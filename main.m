clc;clear;
fundamp = 57.73;%fundamental amplitude*********1
amp_mod = 0;%amplitude modulation depth**********2
fre_mod = 2;%modulation frequency**********3
iniang_mod = 0;%modulation phase angle**********4
deltaf = 0;%frequency offset**********5
ang_mod = 0;%phase modulation depth**********6
dfre = 0;%rate of change of frequency**********7
iniang = pi/2;%initial phase angle**********8
para = [fundamp, amp_mod, fre_mod, iniang_mod, deltaf, ang_mod, dfre, iniang];

fs = 1000;%sampling frequency
tall = 5;%test duration
t = 0: 1/fs: tall-1/fs;
fr = 1000;%PMU reporting rate
dist = fs/fr;
[sig, ref_val] = get_signal(para, t, dist);%generate test signals
%step test
% para1 = para;
% para1(1) = 1.1*fundamp;%amplitude step
% % para1(8) = iniang + 10*pi/180;%phase step
% [sig1, ref_val1] = get_signal(para1, t, dist);
% sig = [sig, sig1];
% ref_val = [ref_val; ref_val1];
% tall = 2*tall;

% t = 0: 1/fs: tall-1/fs;
% fh = (50+deltaf)*2;
% sig = sig + sqrt(2)*0.1*fundamp*cos(2*pi*fh*t + iniang);%add harmonics
% sig = awgn(sig, 70, 'measured');%add noise

fs = 2000;
h1 = load('lpf1.mat'); h1 = h1.h;
winlen = 0.02;
t = -winlen/2: 1/fs: winlen/2;
t = t(2: 2: end);
cp = exp(-1i*2*pi*50*t);
h1 = h1 .* cp;%shift frequency: lowpass to bandpass

h2 = load('lpf2.mat'); h2 = h2.h;
t = (-winlen/2-2/fs): 1/fs: (winlen/2+2/fs);
t = t(2: 2: end);
cp = exp(-1i*2*pi*50*t);
h2 = h2 .* cp;%shift frequency: lowpass to bandpass
hk1 = conv(h1, h2);%the first complex bandpass filter

h3 = load('lpf3.mat'); h3 = h3.h;
cp = exp(-1i*2*pi*60*t);
h3 = h3 .* cp;%shift frequency: lowpass to bandpass
hk2 = conv(h1, h3);%the second complex bandpass filter

fs = 1000;
N = round(winlen*fs);

spec = fft(hk1', fs*100);
mag = abs(spec);  hk1 = hk1 / max(mag); mag1 = mag/max(mag); mag1 = mag1(4801: 5201);
spec = fft(hk2', fs*100);
mag = abs(spec);  hk2 = hk2 / max(mag); mag2 = mag/max(mag); mag2 = mag2(4801: 5201);
mrt = mag2 ./ mag1;
ff = 48: 0.01: 52; %figure; plot(ff, mrt);
p = polyfit(mrt, ff', 2);%polynomial fitting

num = fr * tall;
fre = zeros(num, 1);
dfre = zeros(num, 1);
rt = zeros(num, 1);

idx = 0;
for k=0: num
    cp = k*dist + 1;
    idx = idx + 1;
    if cp-N<=0 || cp+N>length(sig)
        continue;
    end
    t0 = mod(k, fr)/ fr;    
    cal_val = sig(cp-N: cp+N);
    fpha = sum(hk1 .* cal_val);
    amp1 = abs(fpha) * sqrt(2);
    fpha = sum(hk2 .* cal_val);
    amp2 = abs(fpha) * sqrt(2);
    rt(k) = amp2/amp1;
    
    fre(idx) = p(1)*rt(k)^2 + p(2)*rt(k) + p(3);%frequency estimation
end
flen = 1; 
for k=2*flen+1: num
    dfre(k) = (fre(k)-fre(k-2*flen))/(2*flen)*fr;%calculate ROCOF
end

%calculate estimation errors
rfre = ref_val(:, 3); rdf = ref_val(:, 4);
efre = rfre - fre; edf = rdf - dfre;
maxfre = max(abs(efre(fr: num-fr))); maxdf = max(abs(edf(fr: num-fr)));
res = [maxfre, maxdf]; resval = [efre, edf];





