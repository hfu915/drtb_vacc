%% goveqs_scaleup_ct1.m
% construct governing equations for linear scale-up

function out = goveqs_scaleup_ct1(t, in, M0, M1, times, i, s, r, p)

scale = max(min((t-times(1))/(times(2)-times(1)),1),0);
Mt = M1; Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
out = goveqs_basis_ct1(t, in, Mt, i, s, r, p);
