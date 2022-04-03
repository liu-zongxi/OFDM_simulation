%-----------------------生成N个调制信号--------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月31日-----------------%
% 输入
% npsk: 调制阶数
% n: 需要生成的个数
% 输出
% modulated_symbols：生成的N个调制信号,(N,1)
function modulated_symbols = ModSymbolGenerator(npsk, n)
As = [1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];
M = 2^npsk;
A = As(npsk);
data = randi([0, M-1], n, 1);
modulated_symbols = qammod(data, M)./A;
