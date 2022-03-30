%-----------------------LS信道估计----------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月28日-----------------%
function output_H = ChannelEstimation_LS(input_TD, X_pilot, npilot, l_pilot, start_pilot, nused, nframe)
% 输入
% input_TD: (Nused+Npilot, Nframe)，输入的含导频的数据信号
% X_pilot：(Npilot, Nframe) 原始的导频信号
% l_pilot：导频间隔长度
% start_pilot：导频起始位置
% npilot：导频共计多少个这与Nk和l_pilot息息相关
% nframe: 一个发送多少个OFDM符号
% nused: 一个OFDM符号发送了多少数据
% 输出
% output_H: (Nused+Npilot, Nframe),插值后的信道估计结果

output_H = zeros(nused+npilot, nframe);
% 首先，我们要把导频从input_TD中提取出来
indexs_pilot = start_pilot+(0:npilot-1).*l_pilot;
Y_pilot = input_TD(indexs_pilot, :);
% 估计，书中公式(6.7)
H_pilot_LS = Y_pilot ./ X_pilot;
% 下面 该插值了,插值结果是数据+导频的长度
indexs_H = 1 : nused+npilot;
% 此时不可以直接进行插值，因此插值只能插中间，不能插两头，我们要人为的把两头算出来补上
[H_pilot_compensation, indexs_pilot_compensation]=InterpolateCompensation(H_pilot_LS, indexs_pilot, indexs_H(end), nframe);
% 既然补偿的时候就是用的线性插值，就别尝试spline了，都一样
for kk = 1:nframe
    output_H(:, kk) = interp1(indexs_pilot_compensation, H_pilot_compensation(:, kk), indexs_H, "linear");
end
end
