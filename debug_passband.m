%-------------仿真OFDM限幅后的PAPR和BER---------------------%
%-----------------------author:lzx-------------------------%
%-------------------date:2022年4月1日-----------------------%
%% 设置参数
clear;clc;clf;
Nk = 128;           % 子载波个数
Nfft = 128;          % fft长度
CR = 1.2;
Nmod = 2;              % 调制符号所含比特
M = 2^Nmod;         % 调制阶数
NGI = 32;           % 保护间隔长度
fs = 2e6;           % 采样频率
Nos = 8;               % 过采样系数
Tsym = 1/(fs/Nk);   % 系统周期
Ts = 1/(fs*Nos);      % 采样周期
fc = 2e6;           % 载波频率
wc = 2*pi*fc;       % 载波角频率
Nsym = Nfft+NGI;    % 系统长度
t = [0:Ts:2*Tsym-Ts].'/Tsym;
f=[0:fs/(Nfft*2):Nos*fs-fs/(Nfft*2)]-Nos*fs/2;
t0 = t((Nfft/2-NGI)*Nos);   % 把CP开始的时刻作为0时刻
% 滤波器
Fs = 8;             % 滤波器采样频率
Norder = 104;        % 滤波器阶数
dens=20;            % Density factor of filter
FF=[0 1.4 1.5 2.5 2.6 Fs/2]; % Stopband/Passband/Stopband frequency edge vector
WW=[10 1 10]; % Stopband/Passband/Stopband weight vector
h = firpm(Norder,FF/(Fs/2),[0 0 1 1 0 0],WW,{dens});
%%
X_mod = ModSymbolGenerator(Nmod, Nk);       % 生成信号，Nk个
X_mod(1) = 0;                               % 去除直流分量
x = IFFTOversampling(X_mod, Nfft, Nos);     % ifft+过采样，得到的时域信号长度为Nfft*Nos,输入时频域分居两侧
x_GI = AddGI(x, Nfft*Nos, NGI*Nos, 1, "CP");% 添加GI，也要过采样
x_os = [zeros((Nfft/2-NGI)*Nos,1); x_GI; zeros(Nfft*Nos/2,1)];% 这不是过采样，而是为了展示CP而不得已看两个周期的x而做的补零
x_passband = real(x_os.*exp(1j*2*wc*t));    % 上变频
x_clipped = Clipping(x_passband,CR);        % 限幅
X_filter= filter(h,1,fft(x_clipped, Nos*Nfft));       % 先做fft，再做滤波
x_filter = ifft(X_filter);                  % 变回去
x_baseband = x_filter.*exp(-1j*2*wc*t);     % 下变频     

%% 
figure(1); clf % Fig. 7.15(a), (b)
% nn是包含了CP和x的长度,nn1是CP,nn2是信号
nn=(Nfft/2-NGI)*Nos+[1:Nfft*Nos].'; nn1=Nfft/2*Nos+[-NGI*Nos+1:0].'; nn2=Nfft/2*Nos+[0:Nfft*Nos].';
subplot(221)
plot(t(nn1)-t0, abs(x_os(nn1)),'k:'); hold on;
plot(t(nn2)-t0, abs(x_os(nn2)),'k-');
axis([t([nn1(1) nn2(end)])-t0;  0;  max(abs(x_os),[], 1)]);
title(['Baseband signal, with CP']);
xlabel('t (normalized by symbol duration)'); ylabel('abs(x''[m])');
subplot(223)
X_os_dB = 20*log10(abs(fft(x_os)));
plot(f,fftshift(X_os_dB)-max(X_os_dB, [], 1),'k');
xlabel('frequency[Hz]'); ylabel('PSD[dB]'); axis([f([1 end]) -100 0]);
subplot(222)
h_x_passband=histogram(x_passband(nn),50);
xlabel('x'); ylabel('pdf'); title(['Unclipped passband signal']);
subplot(224)
X_passband_dB = 20*log10(abs(fft(x_passband)));
plot(f,fftshift(X_passband_dB)-max(X_passband_dB),'k');
xlabel('frequency[Hz]'); ylabel('PSD[dB]'); axis([f([1 end]) -100 0]);