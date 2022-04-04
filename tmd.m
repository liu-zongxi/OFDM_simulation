%% 设置参数
clear; clf; clc;
% 信噪比
SNRdBs=[10:10];
N_SNR=length(SNRdBs);
% OFDM相关
Maxiter = 1;          % 信噪比迭代次数
CRs=[0.8:0.2:1.6];      % 限幅系数，他会被用于RMS的相乘限幅
N_CR=length(CRs);
gss='*^<sd>v';
Nmod = 2;               % 映射
M = 2^Nmod;             % 调制阶数
Nk = 128;               % 子载波次数，一般Nk=Nfft
Nfft = 128; 
NGI = 32;               % GI长度
Nframe = 1;             % 一帧6个符号
Nsym = Nfft + NGI;      % 系统长度
% 采样
fs = 1e6;               % 采样次数，也就是两个点直接距离
Nos = 8;                % 过采样倍数
Tsym = 1/(fs/Nfft);     % Nsym占了多长时间？
Ts = 1/(fs*Nos);        % 采样周期
ts = [0:Ts:Nframe*Tsym-Ts].';
% 载波
fc = 2e6;               % 上变频载波频率
wc = 2*pi*fc;           % 角频率
% 滤波器
Fs=8; Norder=104; dens=20; 
FF=[0 1.4 1.5 2.5 2.6 Fs/2];
WW=[10 1 10];
h = firpm(Norder,FF/(Fs/2),[0 0 1 1 0 0],WW,{dens});
% 初始化一些储存数组
CF = zeros(Maxiter,1);
% 一些位置
% 生成一帧数据，串并转换，并QPSK，生成一帧
Frame_FDserial = ones(1, Nk*Nframe*Nmod);
% Frame_FDserial = rand(1,Nk*Nframe*Nmod) > 0.5;     % 发送的是bit
Frame_FDparallel = reshape(Frame_FDserial,Nk,Nframe*Nmod);% 串并转换
Frame_mod = QPSKMod(Frame_FDparallel,Nk,Nframe);     %调制
Frame_oversampling = [zeros(1, Nframe); Frame_mod(Nk/2+1:Nk,:); zeros(Nk*(Nos-1)-1,Nframe); Frame_mod(1:Nk/2,:)];
frame = ifft(Frame_oversampling, Nfft*Nos);
frame_serial = reshape(frame, Nframe*Nfft*Nos, 1);
frame_passband =  sqrt(2)*real(frame_serial.*exp(1j*wc*ts));
frame_baseband =  sqrt(2).*frame_passband.*exp(-1j*wc*ts)
frame_parallel = reshape(frame_baseband, Nfft*Nos, Nframe);
Frame_after = fft(frame_parallel, Nfft*Nos)
