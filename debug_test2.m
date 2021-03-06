clear;clc;
input = ones(128,6);
input_TD = ifft(ifftshift(input), 128);
h = [1 2 3 4 5];
input_serial = reshape(input_TD, 1, 128*6);
input_conv = conv(input_serial, h);
input_fading = input_conv(:, 1:length(input_serial));
output_parallel = reshape(input_fading, 128, 6);
output_FD = fftshift(fft(output_parallel, 128));
H = fftshift(fft([h zeros(1, 128-length(h))].', 128));
output_eq = output_FD ./ repmat(H, 1, 6)
n_biterror = sum(sum(output_eq-input~=0))