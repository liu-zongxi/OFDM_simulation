%-----------------------添加GI-----------------------------%
%-----------------------author:lzx-------------------------%
function output_TD = AddGI(input_TD, nfft, nGI, nframe, type_GI)
if type_GI=="CP"    % 实现CP
    output_TD = [input_TD(nfft-nGI+1:nfft, :); input_TD(1:nfft, :)];
elseif type_GI=="ZP" % 实现ZP
    output_TD = [zeros(nGI,nframe); input_TD(1:nfft, :)];
end