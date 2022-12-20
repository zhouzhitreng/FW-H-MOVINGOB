function [t1_new_arr,t1_old_arr,nt1]=pre_intep(t1_new,t)
[t1_new_s,i_s]=min(abs(t1_new));  %这里默认了初始时刻必须为0
i_s=i_s+1;
lent=length(t);
nt1=lent-i_s+1;
t1_old_arr=t;
t1_new_arr=t1_new(end-nt1+1:end);
end

