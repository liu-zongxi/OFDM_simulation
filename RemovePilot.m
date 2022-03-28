%-----------------------在数据中去除导频-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月28日22:26:00--------%
function output_data = RemovePilot(input_data, l_pilot, start_pilot, npilot)

indexs_pilot = start_pilot+(0:npilot-1).*l_pilot;
input_data(indexs_pilot, :) = [];
output_data = input_data;
