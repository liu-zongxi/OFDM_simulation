%-----------------------解调QPSK---------------------------%
%-----------------------author:lzx-------------------------%
function outputs = QPSKDemod(input_data,nk,nframe)
% 输入
% input_data: (Nk, Nframe), 一个频域的复数，会被拆开解调
% nk: 频域并联
% nframe: 一帧包含符号数
% 输出
% outputs：(Nk, 2*Nframe), 解调后，多出一倍，全是01
outputs = zeros(nk,2*nframe);
A = 1/sqrt(2); 
input_data = input_data ./ A;
outputs((1:nk),(1:2:2*nframe-1)) = real(input_data((1:nk),(1:nframe)))>=0;
outputs((1:nk),(2:2:2*nframe)) = imag(input_data((1:nk),(1:nframe)))>=0;