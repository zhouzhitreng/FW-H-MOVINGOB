function [ifft_F1_part,ifft_F2_part,ifft_Q_part]=intep_part(ifft_F1,ifft_F2,ifft_Q,t1_new_arr,t1_old_arr,nt1)
ifft_F1_part=interp1(t1_old_arr,ifft_F1,t1_new_arr,'spline');
ifft_F2_part=interp1(t1_old_arr,ifft_F2,t1_new_arr,'spline');
ifft_Q_part=interp1(t1_old_arr,ifft_Q,t1_new_arr,'spline');

end

