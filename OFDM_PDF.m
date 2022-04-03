%--------------------------OFDM信道的PDF--------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月31日-----------------%

%% 参数设置
clear; clc;
N = 32;          % 生成的调制符号的个数
Npsk = 2;       % 调制的比特
M = 2^Npsk;     % 阶数
Nos = 16;       % 过采样倍数
Nfft = N*Nos;   % fft采样点
Ts = 1/Nfft;    % 采样周期
ts = 0:Ts:1-Ts; % 采样时间
Nhist = 1e3;    % 统计重复次数
Nbin = 30;      % 画直方图时被分成多少类
%% 主程序
for jj = 1:Nhist
    X_mod = ModSymbolGenerator(Npsk, N);
    X(1) = 0+1j*0;
    xI = zeros(Nfft, N);
    xQ = zeros(Nfft, N);
    for ii = 1:N
        % 这是在做什么？计算一下就知道啦
        % 为了作图，每个FFT只使用一个子载波
        % 用那些呢？最左边的和最右边的
        % 为什么？虚拟子载波吧？不太清楚
       if ii<=N/2  
           x = ifft([zeros(1,ii-1) X_mod(ii) zeros(1,Nfft-ii+1)],Nfft);
       else  
           x = ifft([zeros(1,Nfft-N+ii-1) X_mod(ii) zeros(1,N-ii)],Nfft);
       end
       xI(:, ii) = real(x); 
       xQ(:, ii) = imag(x);
    end
    HistI((Nfft*(jj-1)+1):(Nfft*jj)) = sum(xI,2);
    HistQ((Nfft*(jj-1)+1):(Nfft*jj)) = sum(xQ,2);
end
sum_xI = sum(xI,2); sum_xQ = sum(xQ,2); % 这组合出来的就是OFDM符号啦
figure(1), clf, subplot(311);
plot(ts,xI,'k:','linewidth',1), hold on, plot(ts,sum_xI,'b','linewidth',2);
subplot(312)
plot(ts,xQ,'k:','linewidth',1); hold on, plot(ts,sum_xQ,'b','linewidth',2);
subplot(313), plot(ts,abs(sum_xI+1j*sum_xQ),'b','linewidth',2); hold on;
figure(2), clf, subplot(311)
sqrt(var(HistI))
h_xI = histogram(HistI,Nbin);
subplot(312)
h_xQ = histogram(HistQ,Nbin);
subplot(313)
sqrt(var(HistQ))
sqrt(1/(2*N))/Nos
h_x = histogram(abs(HistI+1j*HistI),Nbin);