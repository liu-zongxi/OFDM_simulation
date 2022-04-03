%--------------------对通频带信号限幅------------------------%
%-----------------------author:lzx-------------------------%
function [x_clipped,sigma]=Clipping(x,CL,sigma)
% 输入
% x:要被限幅的信号
% CL: 限幅倍数
% sigma:信号的RMS
% 输出
% x_clipped:限幅后的信号
% sigma:信号的RMS

% 先计算信号的RMS
if nargin<3
  x_mean=mean(x);
  x_dev=x-x_mean;
  sigma=sqrt(x_dev' * x_dev/length(x));
end
% 参考书7.19
A = CL*sigma;
x_clipped = x;
ind = find(abs(x)>A);
x_clipped(ind) = x(ind)./abs(x(ind))*A;