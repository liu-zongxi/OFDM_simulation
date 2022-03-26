%-----------------------去除CP----------------------------%
%-----------------------author:lzx-------------------------%
function output_TD = RemoveGI(input_TD,nfft,nGI)
% 输入
% input_TD: (Nsym,Nframe)输入的并联时域数据
% nfft：fft长度
% nGI: GI长度
% 输出
% output_TD: (Nfft,Nframe)去掉GI后的并联时域数据
    output_TD = input_TD(nGI+1:nfft+nGI,:);
