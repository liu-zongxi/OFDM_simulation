%-----------------------QPSK实现----------------------------%
%-----------------------author:lzx-------------------------%
function outs = QPSKMod(input_data,nk,nframe,m)
% 输入
% input_data: 待数字调制的输入数据(Nk,Nframe*M)
% nk: 子载波个数，也就是并联个数
% nframe: 一帧中包含多少OFDM符号
% m: 调制数
% 输出
% outs: (nk,nframe),输出a+bi
% out_coss:(nk,nframe)，输出实部
% out_sins:(nk,nframe)，输出虚部

if nargin < 4                   % 设置默认值
    m = 2;
end

outs = zeros(nk,nframe);    % 初始化
out_cos = zeros(nk,nframe);
out_sin = zeros(nk,nframe);

input_data_pn = 2*input_data - 1;     % 把0，1映射到-1，1
A = 1/sqrt(2);                        % 归一化幅值

out_cos((1:nk),(1:nframe)) = input_data_pn((1:nk), (1:2:2*nframe-1)) .* A;    % 每次使用两列
out_sin((1:nk),(1:nframe)) = input_data_pn((1:nk), ((2:2:2*nframe))) .* A;

outs = out_cos + out_sin * 1j;