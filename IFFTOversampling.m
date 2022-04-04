%-------------仿真OFDM限幅后的PAPR和BER---------------------%
%-----------------------author:lzx-------------------------%
%-------------------date:2022年4月3日-----------------------%
function [output_TD, time] = IFFTOversampling(input_data, n, nos)
% 输入
% input_data: 输入的信号（N，Nframe）
% n: 原本应该的Nfft=Nk
% nos: 过采样倍数
% 输出
% output_TD:(n*nos,1)输出的时域信号
% time:过采样的时间点

% 如果们没有输入Nos, 就不进行过采样
if nargin<3,  nos=1;  end
nframe = size(input_data, 2);
nfft = n*nos;
Ts = 1/nfft;
time = 0:Ts:1-Ts;
% 为什么要×Nos? 为了保持能量不变
% 这个分居两侧已经包含了ifftshift
output_TD = nos*ifft([input_data(1:n/2, :);  zeros(nfft-n, nframe);  input_data(n/2+1:end, :)], nfft);
