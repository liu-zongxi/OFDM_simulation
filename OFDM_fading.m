%-----------------------用TDL仿真多径衰落-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月25日-----------------%
%% 设置参数
clear;clc;
Nk = 128;           % 子载波个数
Nfft = 128;          % fft长度
Nframe = 6;         % 一帧中有几个OFDM符号
Npsk = 2;              % 调制符号所含比特
M = 2^Npsk;          % 调制数
SR = 250000;        % 符号速率
BR = SR .* Npsk;       % 比特率
NGI = 32;           % 保护间隔长度
EbN0s = 0:2:20;      % 信噪比
Nsym = Nfft+NGI;    % 系统长度
bers = zeros(1,length(EbN0s));  % 误码率储存数组
PowerTDL_dB = [0 -8 -17 -21 -25];   % TDL中信道抽头的功率,dB为单位
Delay = [0 3 5 6 8];                % TDL中信道时延
PowerTDL = 10.^(PowerTDL_dB/10);    % TDL中信道抽头的功率
Nchannel=length(PowerTDL_dB);       % 信道抽头数
Tau_maxTDL = Delay(end)+1;          % 最大时延除以帧长,就是归一化的最大时延
fprintf('EbN0 \t \t ber\t\t\t per\t\t\t nloop \t\t \n');
%% 函数主体

for kk = 1:length(EbN0s)
    % rng('default')          % 初始化随机种子
    EbN0 = EbN0s(kk);
    nloop = 10000;          % 发送多少帧
    n_biterror = 0;         % 错误的数据
    n_bitdata = 0;          % 一共发送了多少数据
    n_packeterror = 0;      % 有多少错误帧
    n_packetdata = 0;       % 发送了多少帧
    for ii = 1:nloop
%--------------------------发射端-------------------------------%
        % 生成一帧数据，串并转换，并QPSK，生成一帧
        frame_FDserial = rand(1,Nk*Nframe*Npsk) > 0.5;     % 发送的是bit
        frame_FDparallel = reshape(frame_FDserial,Nk,Nframe*Npsk);% 串并转换
        frame_mod = QPSKMod(frame_FDparallel,Nk,Nframe);     %调制
        % IFFT
        power_FT = sum(sum(abs(frame_mod).^2))/Nk/Nframe;  % 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_mod);         % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft);             % ifft
        % frame_ifft = ifft(frame_mod, Nfft);
        power_TD = sum(sum(abs(frame_ifft).^2))/Nk/Nframe; % 计算下IFFT前的能量，DT表示时域
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
        power_TDserial = sum(abs(frame_TDserial).^2)/Nk/Nframe;     % 计算出的能量和理论不符啊，没发现问题在哪
        EsN0 = EbN0 + 10*log10(Npsk);                                  % 根据信噪比计算噪声能量，幅值，然后加在信号上
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
        % 信道均衡
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft));
        frame_equalization = frame_recieved_FD ./ repmat(H, 1, Nframe);
        % QPSK解调
 
        frame_demod = QPSKDemod(frame_equalization, Nk, Nframe);
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk*Nframe*Npsk);
       
        % 计算error
        n_biterror_tmp = sum(abs(frame_output-frame_FDserial));
        n_bitdata_tmp = length(frame_FDserial);
        n_biterror = n_biterror + n_biterror_tmp;
        n_bitdata = n_bitdata + n_bitdata_tmp;
        if n_biterror_tmp ~= 0
            n_packeterror = n_packeterror + 1;
        end
        n_packetdata = n_packetdata + 1;
    end
    % 计算在当前信噪比下的误码率
    per = n_packeterror/n_packetdata;
    ber = n_biterror/n_bitdata;
    bers(kk)=ber;
    fprintf('%f\t%e\t%e\t%d\t\n',EbN0,ber,per,nloop);
end
save("BERofdm_rayleigh.mat",'bers');
%% 画图
% bers = load(!"BERofdm_rayleigh.mat").bers;
rayleigh_theory = 0.5.*(1-(1-(1./(10.^(EbN0s./10).*(48/80)+1))).^0.5);
semilogy(EbN0s,rayleigh_theory,'-*',EbN0s,bers,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('理论曲线','实验曲线');