%--------------------------OFDM通频带仿真--------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年4月3日------------------%
%% 设置参数
clear; clf; clc;
% 信噪比
SNRdBs=[0:10];
N_SNR=length(SNRdBs);
% OFDM相关
Maxiter = 1000;          % 信噪比迭代次数
CRs=[0.8:0.2:1.6];      % 限幅系数，他会被用于RMS的相乘限幅
N_CR=length(CRs);
gss='*^<sd>v';
Nmod = 2;               % 映射
M = 2^Nmod;             % 调制阶数
Nk = 128;               % 子载波次数，一般Nk=Nfft
Nfft = 128; 
NGI = 32;               % GI长度
Nframe = 6;             % 一帧6个符号
Nsym = Nfft + NGI;      % 系统长度
% 采样
fs = 1e6;               % 采样次数，也就是两个点直接距离
Nos = 8;                % 过采样倍数
Tsym = 1/(fs/Nsym);     % Nsym占了多长时间？
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
CF = zeros(1, Maxiter);
% 一些位置
bers = zeros(N_CR, N_SNR);
% Z相关
Z_dBs = [2:0.1:16];
N_Z = length(Z_dBs);
CCDF = zeros(1, N_Z);
CCDF_CR = zeros(N_CR, N_Z);
%% 主函数
for ii = 1:N_SNR
    SNR = SNRdBs(ii);
    for kk = 1:N_CR
        CR = CRs(kk);
        nobe = 0;
        for jj = 1:Maxiter
            % 生成一帧数据，串并转换，并QPSK，生成一帧
            Frame_FDserial = rand(1,Nk*Nframe*Nmod) > 0.5;     % 发送的是bit
            Frame_FDparallel = reshape(Frame_FDserial,Nk,Nframe*Nmod);% 串并转换
            Frame_mod = QPSKMod(Frame_FDparallel,Nk,Nframe);     %调制
            % 做过采样ifft
            % https://blog.csdn.net/sinat_38151275/article/details/85268026
            Frame_oversampling = [zeros(1, Nframe); Frame_mod(Nk/2+1:Nk,:); zeros(Nk*(Nos-1)-1,Nframe); Frame_mod(1:Nk/2,:)];
            frame = ifft(Frame_oversampling, Nfft*Nos);
%             Frame_mod(1,:) = 0+1j*0;                                    % 没有直流分量
%             frame = IFFTOversampling(Frame_mod, Nk, Nos);           % 过采样IFFT
            frame_GI = AddGI(frame, Nfft*Nos, NGI*Nos, Nframe, "CP");  % 添加保护间隔
            frame_serial = reshape(frame_GI, Nframe*Nsym*Nos, 1);   % 程序写的有点问题，只能都按列来了
            % frame_padding = [zeros((Nfft/2-NGI)*Nos, Nframe); frame_GI; zeros(Nfft*Nos/2, Nframe)];
            frame_passband = sqrt(2) .* real(frame_serial.*exp(1j*wc*ts));
            % frame_clipped = frame_passband;
            frame_clipped = Clipping(frame_passband, CR);
            % frame_filter = frame_clipped;
            frame_filter = ifft(abs(fft(h.',size(frame_clipped, 1))).*fft(frame_clipped));      % 滤波原理
            if ii == N_SNR
                CF(1,jj) = PAPR_dB(frame_filter);
            end
            % 加噪声
            % frame_noise = frame_filter;
            frame_noise = awgn(frame_filter,SNR,'measured');   % add Noise(AWGN)
            frame_baseband = sqrt(2) .* frame_noise.*exp(-1j*wc*ts);
            % frame_sample = frame_baseband(1+[0:Nos:Nsym*Nos*Nframe-1],1);
            frame_parallel = reshape(frame_baseband, Nsym*Nos, Nframe);
            frame_noGI = RemoveGI(frame_parallel, Nfft*Nos, NGI*Nos);
            Frame_os = fft(frame_noGI,  Nfft*Nos);
            Frame = [Frame_os(Nfft/2+Nfft*(Nos-1)+1:Nfft*Nos, :); Frame_os(2:Nfft/2+1, :)];
            % Frame = fft(frame_parallel, Nfft);
%             for kkk = 1: Nframe
%                 Frame(:,kkk) = fftshift(Frame(:,kkk));
%             end
            Frame_demod = QPSKDemod(Frame, Nfft, Nframe);
            nobe = nobe + biterr(Frame_FDparallel(2:end,:),Frame_demod(2:end,:));

            % nobe = nobe + biterr(Frame_FDparallel(2:end,:),Frame_demod(2:end,:)); 
        end
        bers(kk, ii) = nobe/Maxiter/Nmod/Nk
        if ii == N_SNR
            for jjj = 1:N_Z
                CCDF(jjj) = sum(CF>Z_dBs(jjj))/Maxiter;
            end
            CCDF_CR(kk, :) = CCDF;
        end
    end
end
%% 画图
% figure(1)
% plot(abs(fft(frame_clipped)));
% figure(2)
% plot(abs(fft(frame_filter))); hold off
for kk = 1:N_CR
   gs = gss(kk);
   subplot(211), semilogy(Z_dBs,CCDF_CR(kk,:),[gs '-']), hold on
   subplot(212), semilogy(SNRdBs,bers(kk,:),[gs '-']), hold on
end 