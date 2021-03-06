clear;clc;
Nk = 128;           % 子载波个数
Nfft = 128;          % fft长度
Nframe = 6;         % 一帧中有几个OFDM符号
M = 2;              % 调制符号所含比特
SR = 250000;        % 符号速率
BR = SR .* M;       % 比特率
NGI = 32;           % 保护间隔长度
EbN0s = 10:1:20;      % 信噪比
Nsym = Nfft+NGI;    % 系统长度
PowerTDL_dB = [0 -8 -17 -21 -25];   % TDL中信道抽头的功率,dB为单位
Delay = [0 3 5 6 8];                % TDL中信道时延
PowerTDL = 10.^(PowerTDL_dB/10);    % TDL中信道抽头的功率
Nchannel=length(PowerTDL_dB);       % 信道抽头数
Tau_maxTDL = Delay(end)+1;          % 最大时延除以帧长,就是归一化的最大时延
 % 生成一帧数据，串并转换，并QPSK，生成一帧
        frame_FDserial = rand(1,Nk*Nframe*M) > 0.5;     % 发送的是bit
        frame_FDparallel = reshape(frame_FDserial,Nk,Nframe*M);% 串并转换
        frame_mod = QPSKMod(frame_FDparallel,Nk,Nframe);     %调制
        % IFFT
        power_FT = sum(sum(abs(frame_mod).^2))/Nk/Nframe;  % 计算下IFFT前的能量，FT表示频域
%         frame_mod_shift = ifftshift(frame_mod);         % 频域归零
%         frame_ifft = ifftshift(ifft(frame_mod_shift, Nfft));
        frame_mod_shift = zeros(Nk, Nframe);
        for jj = 1:Nframe
            frame_mod_shift(:,jj) = ifftshift(frame_mod(:,jj));         % 频域归零
        end
        frame_ifft = zeros(Nk, Nframe);
        for jj = 1:Nframe
            frame_ifft = ifft(frame_mod_shift, Nfft);
        end
        % 添加保护间隔
        frame_withGI = AddGI(frame_ifft, Nfft, NGI, Nframe, "CP");  % 添加保护间隔
        % 并串转换
        frame_TDserial = reshape(frame_withGI,1,Nsym*Nframe);
        % 信道衰落
        channel = Rayleigh_model(Nchannel, PowerTDL);
        h = zeros(1, Tau_maxTDL);
        h(Delay+1) = channel;
        frame_conv = conv(frame_TDserial, h);
        frame_fading = frame_conv(:,1:length(frame_TDserial));
        % 接收端，串并转换
        frame_recieved_parallel = reshape(frame_fading,Nsym,Nframe);
        % 去GI
        frame_noGI = RemoveGI(frame_recieved_parallel, Nfft, NGI);
        % FFT
%         frame_noGI_shift = fftshift(frame_ifft);
%         frame_recieved_FD = fftshift(fft(frame_noGI_shift, Nfft));
        frame_recieved_FD = zeros(Nk, Nframe);
        for jj = 1:Nframe
            frame_recieved_FD(:, jj) = fftshift(fft(frame_noGI(:, jj), Nfft));
        end
        % 信道均衡
        H = fftshift(fft([h zeros(1, Nfft-Tau_maxTDL)].', Nfft));
        frame_equalization = frame_recieved_FD ./ repmat(H, 1, Nframe);
        % QPSK解调
        frame_demod = QPSKDemod(frame_equalization, Nk, Nframe);
        % 并串转换
        frame_output = reshape(frame_demod, 1, Nk*Nframe*M);
       
        % 计算error
        n_biterror = sum(abs(frame_output-frame_FDserial))