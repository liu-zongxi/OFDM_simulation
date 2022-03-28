%% 设置参数
clear;clc;
% OFDM相关参数
Nk = 128;           % 子载波个数
Nused = 96;         % 数据个数，剩下来的留给导频了
Nfft = 128;          % fft长度
Nframe = 6;         % 一帧中有几个OFDM符号
Npsk = 2;              % 调制符号所含比特
M = 2^Npsk;          % 调制数
SR = 250000;        % 符号速率
BR = SR .* Npsk;       % 比特率
NGI = 32;           % 保护间隔长度
Nsym = Nfft+NGI;    % 系统长度
% 信噪比相关参数
EbN0s = 0:2:20;      % 信噪比
bers = zeros(1,length(EbN0s));      % 误码率储存数组
% 信道相关参数
PowerTDL_dB = [0 -8 -17 -21 -25];   % TDL中信道抽头的功率,dB为单位
Delay = [0 3 5 6 8];                % TDL中信道时延
PowerTDL = 10.^(PowerTDL_dB/10);    % TDL中信道抽头的功率
Nchannel=length(PowerTDL_dB);       % 信道抽头数
Tau_maxTDL = Delay(end)+1;          % 最大时延除以帧长,就是归一化的最大时延
% 导频信息
L_pilot = 4;            % 导频间隔
start_pilot = 1;        % 导频起始位置
Npilot = Nk/L_pilot;    % 导频数量

%--------------------------发射端-------------------------------%
        % 生成一帧数据，串并转换，并QPSK，生成一帧
        frame_FDserial = rand(1,Nused*Nframe*Npsk) > 0.5;     % 发送的是bit
        frame_FDparallel = reshape(frame_FDserial,Nused,Nframe*Npsk);% 串并转换
        frame_mod = QPSKMod(frame_FDparallel,Nused,Nframe);     %调制
        % 插入导频,梳状导频频率上离散，时间上延续
        [X_pilot, frame_with_pilot] = AddPilot(frame_mod, L_pilot, start_pilot, 1, Npilot, Nused, Nframe);
        % IFFT
        % power_FT = sum(sum(abs(frame_mod).^2))/Nk/Nframe;  % 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_with_pilot);         % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft);             % ifft
        % frame_ifft = ifft(frame_mod, Nfft);
        % power_TD = sum(sum(abs(frame_ifft).^2))/Nk/Nframe; % 计算下IFFT前的能量，DT表示时域
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP");  % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI,1,Nsym*Nframe);
            % x=1:1:160;
            % hold on;
            % plot(x, frame_TDserial(1:160),'b');
%--------------------------Channel-------------------------------%
        % 信号先经历衰落
        channel = Rayleigh_model(Nchannel, PowerTDL);
        h = zeros(1, Tau_maxTDL);
        h(Delay+1) = channel;
        frame_conv = conv(frame_TDserial, h);
        frame_fading = frame_conv(:, 1:length(frame_TDserial));        % 看似是线性卷积，实际上由于CP变成了循环卷积
        % 陈老湿方法添加高斯白噪声，本质上是一样的
%       attn = sqrt(0.5*power_TDserial*SR/BR*10.^(-EbN0/10));
%       noise_msg = attn .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
        frame_recieved = frame_fading;
            % plot(x, noise_msg(1:160),'r');
            % hold off;
%--------------------------接收端-------------------------------%
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved,Nsym,Nframe);
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        % frame_recieved_FD = fft(frame_noGI, Nfft);
        % 信道估计，包含两个部分，导频估计和插值
        H_LS_linear = ChannelEstimation_LS(frame_recieved_FD, X_pilot, Npilot, L_pilot, start_pilot, Nused, Nframe)
        % 信道均衡
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft))
        frame_equalization = frame_recieved_FD ./ repmat(H, 1, Nframe);
        % QPSK解调
 
        frame_demod = QPSKDemod(frame_equalization, Nused, Nframe);
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nused*Nframe*Npsk);
       
        % 计算error
        n_biterror_tmp = sum(abs(frame_output-frame_FDserial))