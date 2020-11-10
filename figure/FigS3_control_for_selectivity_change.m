%% In the name of Allah
%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
%  Fig S3 Object coding in IT

%% Set Paths
clear all
close all
clc

%% Settings
ind_b=[1 500];
Ns=92;
win_=300;step=10;
Percenti=33.33;
control=0;
Shuffle_highlow_PPL=0;
Shuffle_Preferance=1;
baseline_corrected_plv=1;
N_shuff=1001;

inf_cal=@(x,y,z) ROC_(x,y); method='roc';inf_pre=@(x,y) (x-y)./(y+x);

 %% Load
 load ([cd '\Mat\Firing Rate\Normalized_firing_rate.mat']);
load ([cd '\Mat\Firing Rate\firing_rate.mat']);

psth_cr_it_n=Val_cr_IT_psth_norm;psth_wr_it_n=Val_wr_IT_psth_norm;
psth_cr_fef_n=Val_cr_FEF_psth_norm;psth_wr_fef_n=Val_wr_FEF_psth_norm;

psth_cr_it=Val_cr_IT_psth;psth_wr_it=Val_wr_IT_psth;
psth_cr_fef=Val_cr_FEF_psth;psth_wr_fef=Val_wr_FEF_psth;

t_h=[1:2100]-500;
ind_ses=session_selection(1);
idx_it=find(ismember(ind_ses_it,ind_ses));
ind_ses_m1=setdiff(ind_ses,[52:92]);
ind_ses_m2=setdiff(ind_ses,[1:51;]);
it_m1=find(ismember(ind_ses_it,ind_ses_m1));
it_m2=find(ismember(ind_ses_it,ind_ses_m2));
ind_it_m1=(ind_ses_it<52);
ind_it_m2=(ind_ses_it>51);


%
th_h=0.002;p_val_d=0.05;
ixx_it=idx_it;
ih_s=ixx_it(find(abs(mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find(abs(mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ind_it_m1=ind_it_m1(ixx_it);
ind_it_m2=ind_it_m2(ixx_it);

%%
load(strcat(cd,'\Mat\PPL\ROC_test_selectivity.mat'))
x_h=mi_all_out(ixx_it,1);y_h=mi_all_out(ixx_it,2);


figure('name','Fig. S3, Rezayat.E, et al')

hold on
scatter(x_h(ind_it_m1),y_h(ind_it_m1),'b*')
scatter(x_h(ind_it_m2),y_h(ind_it_m2),'gs')
export_file_to_excel([[ones(sum(ind_it_m1),1); 2*ones(sum(ind_it_m2),1)],x_h,y_h],'figS3',8)

y_h=mi_all_out(:,2);x_h=mi_all_out(:,1);
y = y_h;
x = x_h;
% scatter(x_h,y_h)

[r, m ,b]=regression(x_h,y_h,'one');
tbl = table(x,y);
lm = fitlm(tbl,'linear');
xx=-0.3:0.1:1.5;
[r,p]=correlation(x_h,y_h,1001)
length(x_h)
plot(xx,m*xx+b,'r')
xlabel(' Visual selectivity (a.u)')
ylabel('Delay selectivity (a.u)')
ylim([0 1])
xlim([0 1])
set(gca,'fontsize',14,'fontweight','bold')
text( 0.5,0.15,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.5,0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.5,0.05,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')


