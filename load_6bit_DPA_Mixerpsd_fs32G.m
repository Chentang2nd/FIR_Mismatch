clear all;
%----INL DNL----
% rng=200;% set randon number seed
% derinl=0*(rand(1,256)-0.5);% INL with +/- 4 LSB erroer;
% derinl2=0*(rand(1,63)-0.5);
% a=0;
% for i=1:255;
%    rr(i)=a+derinl(i)+1;
%    a=rr(i);
% end,
% a=0;
% for i=1:63;
%    rr2(i)=a+derinl2(i)+1;
%    a=rr2(i);
% end,
% 
% n0=1;
% n1=62;
% p=polyfit(n0:n1,rr2(n0:n1),1);
% INLslop=polyval(p,n0:n1);
% INL=rr2(n0:n1)-INLslop; %INL
% RR6bit(1:63)=rr2(1:63)-32;
% DNL=rr2(n0+1:n1+1)-rr2(n0:n1)-1;
% RR8bit(1:255)=rr(1:1:255)-32;

% figure(300);clf,
% plot(INL,'b-'); hold on;
% plot(DNL,'r-'),grid on;
% grid on;

%-end of ---INL DNL----

Fc=8e9 ;% RF carrier freq
% noise haping parameters of 3rd order errroe feedback
aa=-3; 
bb=-aa;
cc=-1;
DC=0;
%-------FIR---
% N=3;
% Fs=Fc;% data sampled 
% Fp=200e6;% passband-edge frequency 
% Ap=3; %passband ripple
% Ast=60; %stopband attenuation 
% Rp = (10^(Ap/20) - 1)/(10^(Ap/20) + 1); 
% Rst = 10^(-Ast/20);
% NUM = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
% % NUM =floor(10*NUM)/10,
NUM =[1,1];

% fvtool(NUM,'Fs',Fs);
% ----parameter for simulink model----------
AM=1; % amplitude
% bit=3; % # of DS DAC 3 bit
nbit=5 ; % 6bit nequest dac;
g1=0.95/8; %1/8*(1-4/2^4); %unify input amplitude
fs=Fc*4; % sampling frequnecy

osr=8;% %over sampling ratio
ni=2000; % number of initial poit skiped
nn=2^15; % number of FFT point
n=nn;
RBW1=fs/nn,%fft RBW
fbbw=230e6;
fin1=1*fbbw/8;
fin2=2*fbbw/8;
fin3=3*fbbw/8;
fin4=4*fbbw/8;
fin5=5*fbbw/8;
fin6=6*fbbw/8;
fin7=7*fbbw/8;
fin8=8*fbbw/8;

fin1=fs/nn*floor(fin1/fs*nn);
fin2=fs/nn*floor(fin2/fs*nn);
fin3=fs/nn*floor(fin3/fs*nn);
fin4=fs/nn*floor(fin4/fs*nn);
fin5=fs/nn*floor(fin5/fs*nn);
fin6=fs/nn*floor(fin6/fs*nn);
fin7=fs/nn*floor(fin7/fs*nn);
fin8=fs/nn*floor(fin8/fs*nn);
Tres=1/fs/8,
fsrf=fs*8;
RBW2=fsrf/nn/1e6,%fft RBW
figure(50),clf
% for DD=0:1:6;



figure(50),%clf,
% [pY0,f] = pwelch(y3,nn,nn/8,4*nn,fs); %6 bit with DC SD 
 % perfect 6bit

% dbpY0=10*log10(pY0)-max(10*log10(pY0))-41.3;%unify to band psd 
% DD=-3;
% sim('DPA_6bitDAC_MIXER_fs32G');
%     [pxx1,f] = pwelch(y5, nn, nn/8,4*nn,fsrf);
%     dbpxx1=10*log10(pxx1)-max(10*log10(pxx1))-41.3;
% plot(f,dbpxx1,'b*');hold on;
DD=6;
sim('DPA_6bitDAC_MIXER_fs32G');
    [pxx1,f] = pwelch(y5, nn/8, nn/8/8,4*nn/8,fsrf);
    dbpxx1=10*log10( pxx1)-max(10*log10(pxx1))-41.3;
plot(f,dbpxx1,'r-');hold on;
% for DD=4;
% plot(f,dbpxx1,'r-');hold on;
% end
DD=-6;
sim('DPA_6bitDAC_MIXER_fs32G');
[pxx1,f] = pwelch(y5, nn/8, nn/8/8,4*nn/8,fsrf);
dbpxx1=10*log10(pxx1)-max(10*log10(pxx1))-41.3;
plot(f,dbpxx1,'b-');hold on;

% plot(f,dbpY0,'r*-'); grid on;
% end,
fmask=[-85,-85, -70 -70,-41.3, -41.3, -65,-65,-85,-85];
fmm=[1.6e9,3.79e9, 3.8e9, 5.99e9, 6e9,8.49e9, 8.5e9, 10.59e9,10.6e9,16e9];
plot(fmm, fmask , 'g-');
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (dBm/Mhz)')

% YDC=mean(y5)*31,
% % % end of the .m file  
