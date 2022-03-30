%-----------------------导频和信道均衡----------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月27日21:20:22----------%
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
EbN0s = 1:1:50;      % 信噪比
bers_perfect = zeros(1,length(EbN0s));      % 误码率储存数组
bers_LS = zeros(1,length(EbN0s));      % 误码率储存数组
bers_DFT = zeros(1,length(EbN0s));      % 误码率储存数组
bers_MMSE = zeros(1,length(EbN0s));      % 误码率储存数组
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
fprintf('EbN0 \t \t ber\t\t\t per\t\t\t nloop \t\t \n');
%% 函数主体

for kk = 1:length(EbN0s)
    % rng('default')          % 初始化随机种子
    EbN0 = EbN0s(kk);
    nloop = 10000;          % 发送多少帧
    n_biterror_perfect = 0;         % 错误的数据
    n_biterror_LS = 0;         % 错误的数据
    n_biterror_DFT = 0;         % 错误的数据
    n_biterror_MMSE = 0;         % 错误的数据
    n_bitdata = 0;          % 一共发送了多少数据
    n_packeterror_perfect = 0;      % 有多少错误帧
    n_packeterror_LS = 0;      % 有多少错误帧
    n_packeterror_DFT = 0;      % 有多少错误帧
    n_packeterror_MMSE = 0;      % 有多少错误帧
    n_packetdata = 0;       % 发送了多少帧
    for ii = 1:nloop
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
        frame_fading = frame_conv(:,1:length(frame_TDserial));        % 看似是线性卷积，实际上由于CP变成了循环卷积
        % 添加高斯白噪声
        power_TDserial = sum(abs(frame_TDserial).^2)/Nsym/Nframe;     % 计算出的能量和理论不符啊，没发现问题在哪
        EsN0 = EbN0 + 10*log10(Npsk) + 10*log10(Nused/Nsym);                                  % 根据信噪比计算噪声能量，幅值，然后加在信号上
        N0 = power_TDserial .* 10.^(-EsN0/10);
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
        frame_recieved = frame_fading + noise_msg;
        % 陈老湿方法添加高斯白噪声，本质上是一样的
%       attn = sqrt(0.5*power_TDserial*SR/BR*10.^(-EbN0/10));
%       noise_msg = attn .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
%       frame_recieved = frame_TDserial + noise_msg;
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
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft));
        Rhh = H*H';
        H_LS_linear = ChannelEstimation_LS(frame_recieved_FD, X_pilot, Npilot, L_pilot, start_pilot, Nused, Nframe);
        H_MMSE = ChannelEstimation_MMSE(H_LS_linear, Rhh, Nused, Npilot, Nframe, EbN0);
        % 做DFT滤除噪声
        h_LS = ifft(ifftshift(H_LS_linear),Nfft);
        h_LS_DFT = h_LS(1:Tau_maxTDL, :);
        H_LS_DFT = fftshift(fft(h_LS_DFT, Nfft));
        % 信道均衡
        frame_equalization_perfect = frame_recieved_FD ./ repmat(H, 1, Nframe);
        frame_equalization_LS = frame_recieved_FD ./ H_LS_linear;
        frame_equalization_DFT = frame_recieved_FD ./ H_LS_DFT;
        frame_equalization_MMSE = frame_recieved_FD ./ H_MMSE;
        % 去除导频
        frame_no_pilot_perfect = RemovePilot(frame_equalization_perfect, L_pilot, start_pilot, Npilot);
        frame_no_pilot_LS = RemovePilot(frame_equalization_LS, L_pilot, start_pilot, Npilot);
        frame_no_pilot_DFT = RemovePilot(frame_equalization_DFT, L_pilot, start_pilot, Npilot);
        frame_no_pilot_MMSE = RemovePilot(frame_equalization_MMSE, L_pilot, start_pilot, Npilot);
        % QPSK解调
        frame_demod_perfect = QPSKDemod(frame_no_pilot_perfect, Nused, Nframe);
        frame_demod_LS = QPSKDemod(frame_no_pilot_LS, Nused, Nframe);
        frame_demod_DFT = QPSKDemod(frame_no_pilot_DFT, Nused, Nframe);
        frame_demod_MMSE = QPSKDemod(frame_no_pilot_MMSE, Nused, Nframe);
        % 并串转换
        frame_output_perfect = reshape(frame_demod_perfect, 1, Nused*Nframe*Npsk);
        frame_output_LS = reshape(frame_demod_LS, 1, Nused*Nframe*Npsk);
        frame_output_DFT = reshape(frame_demod_DFT, 1, Nused*Nframe*Npsk);
        frame_output_MMSE = reshape(frame_demod_MMSE, 1, Nused*Nframe*Npsk);
        % 计算error
        n_biterror_perfect_tmp = sum(abs(frame_output_perfect-frame_FDserial));
        n_biterror_LS_tmp = sum(abs(frame_output_LS-frame_FDserial));
        n_biterror_DFT_tmp = sum(abs(frame_output_DFT-frame_FDserial));
        n_biterror_MMSE_tmp = sum(abs(frame_output_MMSE-frame_FDserial));
        n_bitdata_tmp = length(frame_FDserial);
        n_biterror_perfect = n_biterror_perfect + n_biterror_perfect_tmp;
        n_biterror_LS = n_biterror_LS + n_biterror_LS_tmp;
        n_biterror_DFT = n_biterror_DFT + n_biterror_DFT_tmp;
        n_biterror_MMSE = n_biterror_MMSE + n_biterror_MMSE_tmp;
        n_bitdata = n_bitdata + n_bitdata_tmp;
        if n_biterror_perfect_tmp ~= 0
            n_packeterror_perfect = n_packeterror_perfect + 1;
        end
        if n_biterror_LS_tmp ~= 0
            n_packeterror_LS = n_packeterror_LS + 1;
        end
        if n_biterror_DFT_tmp ~= 0
            n_packeterror_DFT = n_packeterror_DFT + 1;
        end
        if n_biterror_MMSE_tmp ~= 0
            n_packeterror_MMSE = n_packeterror_MMSE + 1;
        end
        n_packetdata = n_packetdata + 1;
    end
    % 计算在当前信噪比下的误码率
    per_perfect = n_packeterror_perfect/n_packetdata;
    per_LS = n_packeterror_LS/n_packetdata;
    per_DFT = n_packeterror_DFT/n_packetdata;
    per_MMSE = n_packeterror_MMSE/n_packetdata;
    ber_perfect = n_biterror_perfect/n_bitdata;
    ber_LS = n_biterror_LS/n_bitdata;
    ber_DFT = n_biterror_DFT/n_bitdata;
    ber_MMSE = n_biterror_MMSE/n_bitdata;
    bers_perfect(kk)=ber_perfect;
    bers_LS(kk)=ber_LS;
    bers_DFT(kk)=ber_DFT;
    bers_MMSE(kk)=ber_MMSE;
    fprintf('%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t\n',EbN0,ber_perfect,ber_LS,ber_DFT, ber_MMSE, per_perfect,per_LS,per_DFT, per_MMSE, nloop);
end
save("BERofdm_perfect.mat",'bers_perfect');
save("BERofdm_LS.mat",'bers_LS');
save("BERofdm_DFT.mat",'bers_DFT');
save("BERofdm_MMSE.mat",'bers_MMSE');
%% 画图
% bers = load(!"BERofdm_rayleigh.mat").bers;
EbN0s = 1:1:50;      % 信噪比
bers_perfect = load("BERofdm_perfect.mat").bers_perfect;
bers_LS = load("BERofdm_LS.mat").bers_LS;
bers_DFT = load("BERofdm_DFT.mat").bers_DFT;
bers_MMSE = load("BERofdm_MMSE.mat").bers_MMSE;
rayleigh_theory = 0.5.*(1-(1-(1./(10.^(EbN0s./10).*(96/128)+1))).^0.5);
semilogy(EbN0s,rayleigh_theory,'-*',EbN0s,bers_perfect,'-+',EbN0s,bers_LS,'-^', EbN0s,bers_DFT,'->', EbN0s,bers_MMSE,'-<');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('理论曲线','perfect', "LS", "LS-DFT", "MMSE");