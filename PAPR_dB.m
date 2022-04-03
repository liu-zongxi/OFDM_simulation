%-----------------------计算信号的PAPR,单位dB----------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月31日-----------------%
function [PAPR_dB, power_avg_dB, power_peak_dB] = PAPR_dB(x)
Nx = length(x);
xI = real(x);
xQ = imag(x);

power = xI.*xI + xQ.*xQ;
power_peak = max(power);
power_peak_dB = 10*log10(power_peak);
power_avg = sum(power)/Nx;
power_avg_dB = 10*log10(power_peak/power_avg);
PAPR_dB = 10*log10(power_peak/power_avg);