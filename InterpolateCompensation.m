%------------------补偿导频估计的头尾，便于插值---------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月28日-----------------%
function [output_pilot, indexs_pilot_compensation] = InterpolateCompensation(input_pilot, indexs_pilot, length_H, nframe)
% 输入
% input_pilot: (Nused, Nframe)，输入导频信号
% indexs_pilot：导频所在的位置
% length_H：需要补偿的长度，也就是Nused+Npilot
% nframe: 一个发送多少个OFDM符号

% 输出
% output_pilot：(Npilot+1, Nframe)补偿了头或者尾的导频
% indexs_pilot_compensation： 补偿了头尾的导频位置
head_compensation = zeros(1, nframe);
tail_compensation = zeros(1, nframe);
for kk = 1:nframe
    % 如果前面有空缺，先补偿头，就是第一个数
    if indexs_pilot(1)>1
        slope = (input_pilot(2, kk)-input_pilot(1, kk))/(indexs_pilot(2)-indexs_pilot(1));
        head_compensation(:, kk) = input_pilot(1, kk)-slope*(indexs_pilot(1, kk)-1);
    end
    % 尾部空缺再补偿尾，就是最后一个数
    if  indexs_pilot(end)<length_H
        slope = (input_pilot(end)-input_pilot(end-1))/(indexs_pilot(end)-indexs_pilot(end-1));
        tail_compensation(:, kk) = input_pilot(end)+slope*(length_H-indexs_pilot(end));
    end
end
if indexs_pilot(1)>1
    indexs_pilot_compensation = [1 indexs_pilot];
    output_pilot = [head_compensation; input_pilot];
end
if  indexs_pilot(end)<length_H
    indexs_pilot_compensation= [indexs_pilot length_H];
    output_pilot = [input_pilot; tail_compensation];
end
end