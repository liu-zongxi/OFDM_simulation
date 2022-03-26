%-----------------------生成瑞利信道-----------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月25日-----------------%
function H=Rayleigh_model(nchannel, power_channel)
% 瑞利衰落信道
% 输入
% nchannel： 多径信道的个数
% power_channel：（1, nchannel），每一个信道的功率
% 输出
% H:(1, nchannel),一个瑞利信道,符合高斯分布的nchannel个随机数，代表着衰落
H = (randn(1,nchannel)+1j*randn(1,nchannel)).*sqrt(power_channel/2);
% 功率除以二的原因是瑞利分布的E(x^2)=2\sigma^2