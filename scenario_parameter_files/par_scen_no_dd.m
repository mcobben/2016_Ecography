y_max=100;
x_max=25;
K=100;

N_deme=x_max*y_max;
N_ind=N_deme*K; %250000;

s_migr=0.5;
disp_rate=0.1;
disp_dist=2;
m_rate=10e-5;

repr_factor=0.8;            % reproduction factor at location y_max
rc_rl=(1-repr_factor)/y_max;    % slope of the location correction for reproduction

r_factor=0.4;               % reproduction factor at density=1
rc_rd=1-r_factor;            % slope of the density correction for reproduction

s_factor=1;               % survival factor at density=1
rc_s=1-s_factor;            % slope of the density correction for survival
% perc_std=0.0;               % percentage standard deviation of the local winter survival
std_surv=0;               % std dev of the local winter survival, independent of local winter survival
std_migr=0;             % std dev of the migration decision
std_s_migr=0;         % std dev of the migration survival

