function [x_z,e_z,rt_z] = zscore_by_session(x,e,rt,sess)

x_z = zeros(size(x));
e_z = zeros(size(e));
rt_z = zeros(size(rt));

for ss = unique(sess)'
    idx = sess==ss;
    if sum(idx)>2
        x_z(idx) = zscore(x(idx));
        e_z(idx) = zscore(e(idx));
        rt_z(idx)= zscore(rt(idx));
    end
end
end
