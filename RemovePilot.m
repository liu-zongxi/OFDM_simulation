%-----------------------在数据中去除导频-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月28日22:26:00--------%
function output_data = RemovePilot(input_data, l_pilot, start_pilot, npilot)
% 输入
% input_data: （Nused+Npilot，Nframe）,包含了导频的数据
% l_pilot：导频间隔长度
% start_pilot：导频起始位置
% npilot：导频共计多少个这与Nk和l_pilot息息相关
% nframe: 一个发送多少个OFDM符号
% nused: 一个OFDM符号发送了多少数据
% 输出
% output_data： 不包含导频的纯数据
indexs_pilot = start_pilot+(0:npilot-1).*l_pilot;
input_data(indexs_pilot, :) = [];
output_data = input_data;
