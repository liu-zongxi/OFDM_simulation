%-----------------------OFDM信号的PAPR的CCFD----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月31日-----------------%

%%
clear; clc;
Nffts = 2.^[6:10];          % IFFT点数，这是本仿真的变量
Npsk = 2;                   % 调制对应比特数
M= 2^Npsk;                  % 调制阶数
Nblk = 1e3;                 % 迭代次数
Z_dBs = 4:0.1:10;           % IFFT后时域信号的能量大小
N_Z_dBs = length(Z_dBs);    % 
CCDF_formula = @(N, sigma2, z) 1-((1-exp(-z.^2/(2*sigma2))).^N);    % CCFD计算公式

%% 主函数
for kk = 1:length(Nffts)      % 对于每一个不同的Nfft
    Nfft = Nffts(kk);
    x = zeros(Nfft, Nblk);
    x_CF = zeros(1, Nblk);
    CCDF_simulated = zeros(1,N_Z_dBs);
    for ii = 1:Nblk
        X_mod = ModSymbolGenerator(Npsk, Nfft);
        x(:, ii) = ifft(ifftshift(X_mod), Nfft) .*sqrt(Nfft);
        x_CF(ii) = PAPR_dB(x(:,ii));
    end
    sigma2 = mean(mean(abs(x)))^2/(pi/2);   % E(Z^2) = 2sigma^2
    CCDF_theoretical = CCDF_formula(Nfft, sigma2, 10.^(Z_dBs/20));
    for jj = 1:N_Z_dBs
        CCDF_simulated(jj) = sum(x_CF>Z_dBs(jj))/Nblk;
    end
    semilogy(Z_dBs,CCDF_theoretical,'-');  hold on; grid on;
    semilogy(Z_dBs,CCDF_simulated,':*');
end 
title('OFDM system with N-point FFT');
xlabel('PAPR0[dB]');
ylabel('CCDF=Probability(PAPR>PAPR0)'); 
legend('Theoretical','Simulated');