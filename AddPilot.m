%-----------------------在数据中插入导频-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月27日22:16:56--------%
function [X_pilot, outputs] = AddPilot(input_data, l_pilot, start_pilot, value_pilot, npilot,nused, nframe)
% 输入
% input_data: 
% l_pilot：导频间隔长度
% start_pilot：导频起始位置
% value_pilot：导频是啥，本例中为1
% npilot：导频共计多少个这与Nk和l_pilot息息相关
% nframe: 一个发送多少个OFDM符号
% nused: 一个OFDM符号发送了多少数据
% 输出
% outputs: 被插入了导频的数据
% X_pilot: 导频长啥样，用于信道估计


% 初始化一下outputs
outputs = zeros(nused+npilot, nframe);
% 首先，计算导频要放置的位置
indexs_pilot = start_pilot+(0:npilot-1).*l_pilot;
% 再计算一下数据要放置的位置
indexs_data = 1:nused+npilot;         % 先生成一个长度是数据+导频的数组
indexs_data(indexs_pilot) = [];       % 删掉导频所在位置，就是数据所在位置
% 先给导频赋值
outputs(indexs_pilot, 1:nframe) = value_pilot;
% 再把数据放到该放的位置
outputs(indexs_data, 1:nframe) = input_data;
% 最后，返回一下X_pilot,方便最后的估计， 由于用于估计，他应该添加fftshift
X_pilot = zeros(length(indexs_pilot),nframe);
X_pilot(:,:) = value_pilot;
end