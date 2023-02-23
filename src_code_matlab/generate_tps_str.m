function tps_str = generate_tps_str(tps_int)
% [1, 2, 10] ==> ['001', '002', '010']

tps_str = strings(1, numel(tps_int));
for i = 1:numel(tps_int)
     tps_str(i) = sprintf( '%05d', tps_int(i));
end

end