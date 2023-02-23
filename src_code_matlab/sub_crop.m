function out_data = sub_crop(data, scale)
    out_data = data(scale(1,1):scale(1,2),scale(2,1):scale(2,2),scale(3,1):scale(3,2));
end