%-----------------------MMSE信道估计----------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月29日----------------%
function H_MMSE = ChannelEstimation_MMSE(H_LS, Rhh, nused, npilot, nframe, EbN0)
% 输入
% H_LS: (Nk, Nframe),LS估计的结果,和加权数组相乘就是结果
% Rhh: H自相关矩阵，H*H'
% nframe: 一个发送多少个OFDM符号
% nused: 一个OFDM符号发送了多少数据
% npilot: 导频数量
% 输出
% H_MMSE: 估计结果

% 具体如何计算的请参考https://blog.csdn.net/qq_37989552/article/details/102946707
SNR = 10^(EbN0/10);
H_MMSE = zeros(nused+npilot, nframe);
for kk = 1:nframe
    % Rhh = H_LS(:,kk)*(H_LS(:,kk)');
    W = Rhh/(Rhh+(1/SNR)*eye(nused+npilot));
    % W = Rhh*inv(Rhh+(1/(SNR)).*eye(nused+npilot, nused+npilot));
    H_MMSE(:, kk) = W*H_LS(:,kk);
end