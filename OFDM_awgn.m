%-----------------------仿真OFDM---------------------------%
%-----------------------author:lzx-------------------------%
%% 设置参数
clear;clc;
Nk = 128;           % 子载波个数
Nfft = 128;          % fft长度
Nframe = 6;         % 一帧中有几个OFDM符号
M = 2;              % 调制符号所含比特
SR = 250000;        % 符号速率
BR = SR .* M;       % 比特率
NGI = 32;           % 保护间隔长度
EbN0s = 3:1:10;      % 信噪比
Nsym = Nfft+NGI;    % 系统长度
bers = zeros(1,length(EbN0s));
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
        % 生成一帧数据，串并转换，并QPSK，生成一帧
        frame_FDserial = rand(1,Nk*Nframe*M) > 0.5;     % 发送的是bit
        frame_FDparallel = reshape(frame_FDserial,Nk,Nframe*M);% 串并转换
        frame_mod = QPSKMod(frame_FDparallel,Nk,Nframe);     %调制
        % IFFT
        power_FT = sum(sum(abs(frame_mod).^2))/Nk/Nframe;  % 计算下IFFT前的能量，FT表示频域
        frame_mod_shift = ifftshift(frame_mod);         % 频域归零
        frame_ifft = ifft(frame_mod_shift, Nfft);             % ifft
        power_TD = sum(sum(abs(frame_ifft).^2))/Nk/Nframe; % 计算下IFFT前的能量，DT表示时域
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP");  % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI,1,Nsym*Nframe);
        x=1:1:160;
        % hold on;
        % plot(x, frame_TDserial(1:160),'b');
        % Channel
        power_TDserial = sum(abs(frame_TDserial).^2)/Nsym/Nframe;     % 计算出的能量和理论不符啊，没发现问题在哪
        EsN0 = EbN0 + 10*log10(M);                                  % 根据信噪比计算噪声能量，幅值，然后加在信号上
        N0 = power_TDserial .* 10.^(-EsN0/10);
        noise_msg = sqrt(N0 / 2) .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
        frame_recieved = frame_TDserial + noise_msg;
%         attn = sqrt(0.5*power_TDserial*SR/BR*10.^(-EbN0/10));
%         noise_msg = attn .* (randn(size(frame_TDserial)) + 1j * randn(size(frame_TDserial)));
%         frame_recieved = frame_TDserial + noise_msg;
        % plot(x, noise_msg(1:160),'r');
        % hold off;
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_recieved,Nsym,Nframe);
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
        frame_recieved_FD_shift = fft(frame_noGI, Nfft);
        frame_recieved_FD = fftshift(frame_recieved_FD_shift);
        % QPSK解调
        frame_demod = QPSKDemod(frame_recieved_FD, Nk, Nframe);
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk*Nframe*M);
       
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
save("BERofdm.mat",'bers');
%% 画图
% bers = load("BERofdm.mat").bers;
awgn_theory = [0.0228784075610853,0.0125008180407376,0.00595386714777866,0.00238829078093281,0.000772674815378444,0.000190907774075993,3.36272284196176e-05,3.87210821552205e-06];
semilogy(EbN0s,awgn_theory,'-*',EbN0s,bers,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('理论曲线','实验曲线');