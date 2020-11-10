%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 3
clear all
close all
clc

%% Settings
obj_label = {'Pref' 'Intermediate' 'NPref'};                                 % objects label
loc_label = {'IN' 'OUT'};                                                    % locations label
bands_h_l = { ' \theta' ' \alpha' ' \beta'  '  \gamma' '  H\gamma' '  HH\gamma'};
bands_h_name = { 'theta' 'alpha' 'beta'  'gamma' '  Hgamma' '  HHgamma'};

Color_loc=[1,0,0;0,0,1];
Color_obj=[0,0,0;0.5,0,0.5;0,1,0];
Color_bhv=[0.7,0.3,0;0,0.7,0.7];
freq_bands = [3.9,7.9,21.9,35.9,50,70;8.1,15.1,28.1,50.1,130.1,130.1];


t1_delay= 600; t2_delay= 1200;                                                % Delay period
t1_vis= 50; t2_vis= 350;                                                      % Visual period
t1_base= -300; t2_base= 50;                                                     % Fixation period
t = 1:2100;
t = decimate(t,5,'fir');
t_h = t-500;

%% 
addpath(genpath(strcat(cd,'\Tools\circstat-matlab-master')))
load(strcat(cd,'\Mat\phase_dist\phase_dist_full.mat'))
phase_cr=plv_cr;phase_wr=plv_wr;
load(strcat(cd,'\Mat\PPL\PPL_shuff_1001_2.mat'));


%% prefferance change by multiunit
pref_change =1;ind_nonslective = [];
if pref_change
    load multiunit_preferance_pval.mat
    %     load('prefereance_power.mat')
    %     ind_ses_p = [10 15 18 36 47 50 57 59 85 65 66 72 80 84 67 81 22 43];
    ind_ses_p = [10 15 18 36 47 50 57 80 81 43];
    %     ind_ses_p=[10 15];
    p=0.05;
    ind_nonslective=find((p_rate>p)&(p_baseline>p));
    p2=0.05;
    ind_pow=[];
    %      ind_pow=ind_nonslective(p_val_lfp_power_pref(ind_nonslective)<=p2);
    %      ind_nonslective= ind_nonslective(p_val_lfp_power_pref(ind_nonslective)>p2);
    %      ind_ses_p=[ind_ses_p ind_pow];
    changed_ses=[];
    changed_ses([10],:)=[2 1 3];
    changed_ses([ 15 ],:)=[2 1 3];
    changed_ses([ 81 ],:)=[2 1 3];
    changed_ses([ 43 ],:)=[2 1 3];
    changed_ses([ 66 ],:)=[2 1 3];
    changed_ses([ 80 ],:)=[2 1 3];
    changed_ses([ 84 ],:)=[2 1 3];
    changed_ses([ 72 ],:)=[2 1 3];
    
    changed_ses([ 47 ],:)=[2 1 3];
    changed_ses([ 57 ],:)=[2 1 3];
    changed_ses([59],:)=[2 1 3];
    
    changed_ses([18 ],:)=[1 3 2];
    changed_ses([ 36 ],:)=[1 3 2];
    changed_ses([ 50],:)=[1 3 2];
    
    %     changed_ses(ind_pow,:)=obj_ind_final(ind_pow ,:);
    %         changed_ses([ 92],:)=[1 2 3];
    
    load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Codes\obj_ind_h1.mat')
    
    
    plv_cr_=ppl_cr_n;plv_wr_=ppl_wr_n;
    phase_cr_=phase_cr;phase_wr_=phase_wr;
    
    nochanged_ses=[];
    p=0.01;
    for ss = 1:92
        obj_ix_f=obj_ind_h(ss,:);
        
        %         obj_ix_f=changed_ses(ss,:);
        
        %         plv_cr_(ss,1,:,:) = squeeze(plv_cr(ss,obj_ix_f(1),:,:));
        %         plv_cr_(ss,2,:,:) = squeeze(plv_cr(ss,obj_ix_f(2),:,:));
        %         plv_cr_(ss,3,:,:) = squeeze(plv_cr(ss,obj_ix_f(3),:,:));
        %         plv_cr_(ss,4,:,:) = squeeze(plv_cr(ss,obj_ix_f(1)+3,:,:));
        %         plv_cr_(ss,5,:,:) = squeeze(plv_cr(ss,obj_ix_f(2)+3,:,:));
        %         plv_cr_(ss,6,:,:) = squeeze(plv_cr(ss,obj_ix_f(3)+3,:,:));
        %
        %
        %         plv_wr_(ss,1,:,:) = squeeze(plv_wr(ss,obj_ix_f(1),:,:));
        %         plv_wr_(ss,2,:,:) = squeeze(plv_wr(ss,obj_ix_f(2),:,:));
        %         plv_wr_(ss,3,:,:) = squeeze(plv_wr(ss,obj_ix_f(3),:,:));
        %         plv_wr_(ss,4,:,:) = squeeze(plv_wr(ss,obj_ix_f(1)+3,:,:));
        %         plv_wr_(ss,5,:,:) = squeeze(plv_wr(ss,obj_ix_f(2)+3,:,:));
        %         plv_wr_(ss,6,:,:) = squeeze(plv_wr(ss,obj_ix_f(3)+3,:,:));
        %
        
        phase_cr_(ss,1,:,:) = squeeze(phase_cr(ss,obj_ix_f(1),:,:));
        phase_cr_(ss,2,:,:) = squeeze(phase_cr(ss,obj_ix_f(2),:,:));
        phase_cr_(ss,3,:,:) = squeeze(phase_cr(ss,obj_ix_f(3),:,:));
        phase_cr_(ss,4,:,:) = squeeze(phase_cr(ss,obj_ix_f(1)+3,:,:));
        phase_cr_(ss,5,:,:) = squeeze(phase_cr(ss,obj_ix_f(2)+3,:,:));
        phase_cr_(ss,6,:,:) = squeeze(phase_cr(ss,obj_ix_f(3)+3,:,:));
        
        
        phase_wr_(ss,1,:,:) = squeeze(phase_wr(ss,obj_ix_f(1),:,:));
        phase_wr_(ss,2,:,:) = squeeze(phase_wr(ss,obj_ix_f(2),:,:));
        phase_wr_(ss,3,:,:) = squeeze(phase_wr(ss,obj_ix_f(3),:,:));
        phase_wr_(ss,4,:,:) = squeeze(phase_wr(ss,obj_ix_f(1)+3,:,:));
        phase_wr_(ss,5,:,:) = squeeze(phase_wr(ss,obj_ix_f(2)+3,:,:));
        phase_wr_(ss,6,:,:) = squeeze(phase_wr(ss,obj_ix_f(3)+3,:,:));
        
    end
    
    plv_cr = plv_cr_;plv_wr = plv_wr_;
    phase_cr=phase_cr_;phase_wr=phase_wr_;
    
end

%% session selection
status = 2;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=session_selection_PPL(4);
ind_ses_m2=setdiff(ind_ses,[1:51]);
ind_ses_m1=setdiff(ind_ses,[52:92]);
%%
ppl_cr=plv_cr;ppl_wr=plv_wr;
plv_cr=phase_cr;plv_wr=phase_wr;
NS=92;
% Baseline Correction
baseline_corrected = 1;
ind_b = [41 110];
if baseline_corrected
    ppl_cr_h = ppl_cr;ppl_wr_h = ppl_wr;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_wr = (ppl_wr_h-val_b)./val_sd;
end

%%
f1=freq_bands(1,3);  f2=freq_bands(2,3);
%% Preview Correct
similar_caxis=1;bb1=-0.3;bb2=0.3;bb1_=-2; bb2_=2;
similar_caxis=0;bb1=-0.03;bb2=0.05;bb1_=-0.2; bb2_=0.2;

f=freqs2use;bands_h_l={ ' \theta' ' \alpha' ' \beta' '  \gamma' };

bb1_=-10;bb2_=10;

bins_num=30;

%% average over time and frequency
phase_val_base_cr=[];phase_val_base_wr=[];
phase_val_vis_cr=[];phase_val_vis_wr=[];
phase_val_del_cr=[];phase_val_del_wr=[];

for ss=1:92
    tresh_cr=squeeze(nanmean(nanmean(ppl_cr(ss,1,f>f1&f<f2,:),3),4));
    tresh_wr= squeeze(nanmean(nanmean(ppl_wr(ss,1,f>f1&f<f2,:),3),4));
    tresh_cr=-10;%nanmean([tresh_cr tresh_wr]);
    tresh_wr=tresh_cr;
    ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_cr(ss,1,f>f1&f<f2,t1_base<t_h&t_h<t2_base));
    ppl_h=squeeze(ppl_cr(ss,1,f>f1&f<f2,t1_base<t_h&t_h<t2_base))-squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_base<t_h&t_h<t2_base));
    ind_h=(ppl_h>tresh_cr);
    phase_val_base_cr(ss)=circ_mean( phase_h(ind_h));
    
    
    %  ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_wr(ss,1,f>f1&f<f2,t1_base<t_h&t_h<t2_base));
    %  ppl_h=squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_base<t_h&t_h<t2_base));
    ind_h=(ppl_h>tresh_wr);
    phase_val_base_wr(ss)=circ_mean( phase_h(ind_h));
    
    
    ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_cr(ss,1,f>f1&f<f2,t1_vis<t_h&t_h<t2_vis));
    ppl_h=squeeze(ppl_cr(ss,1,f>f1&f<f2,t1_vis<t_h&t_h<t2_vis))-squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_vis<t_h&t_h<t2_vis));
    ind_h=(ppl_h>tresh_cr);
    phase_val_vis_cr(ss)=circ_mean( phase_h(ind_h));
    
    
    % ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_wr(ss,1,f>f1&f<f2,t1_vis<t_h&t_h<t2_vis));
    %  ppl_h=squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_vis<t_h&t_h<t2_vis));
    ind_h=(ppl_h>tresh_wr);
    phase_val_vis_wr(ss)=circ_mean( phase_h(ind_h));
    
    
    ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_cr(ss,1,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay));
    ppl_h=squeeze(ppl_cr(ss,1,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay))-squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay));
    ind_h=(ppl_h>tresh_cr);
    phase_val_del_cr(ss)=circ_mean( phase_h(ind_h));
    
    
    %  ppl_h=[];phase_h=[];
    phase_h=squeeze(phase_wr(ss,1,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay));
    % ppl_h=squeeze(ppl_wr(ss,1,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay));
    ind_h=(ppl_h>tresh_wr);
    phase_val_del_wr(ss)=circ_mean( phase_h(ind_h));
    
    
    
    
end

%% Fig. 3a Beta band phase differences between areas reflected performance
f=freqs2use;bands_h_l={'\theta' ' \alpha' ' \beta' '\gamma'};
figure('name','Fig. 3a, Rezayat.E, et al')
aa=8;


f1=freq_bands(1,3);  f2=freq_bands(2,3);


ax1=subplot(2,1,1);
%     var_cr=squeeze(plv_cr(ind_ses,1,:,:));
%     var_cr_m1=squeeze(plv_cr(ind_ses_m1,1,:,:));
%     var_cr_m2=squeeze(plv_cr(ind_ses_m2,1,:,:));
%
%      var_ppl_cr1=squeeze(ppl_cr(ind_ses,1,:,:));
%     var_ppl_cr_m1=squeeze(ppl_cr(ind_ses_m1,1,:,:));
%     var_ppl_cr_m2=squeeze(ppl_cr(ind_ses_m2,1,:,:));
%
%     y_m1=(squeeze(circ_mean(circ_mean(var_cr_m1(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
%     y_m2=(squeeze(circ_mean(circ_mean(var_cr_m2(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
%
% %      y_h=(squeeze(circ_mean(circ_mean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),[],3),[],2)));
% %         y_h=(squeeze(circ_median(circ_median(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
%     var_h=var_cr(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay);
%     ppl_h=var_ppl_cr1(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay);
% for ss=1:size(var_h,1)
%     var_h_h=squeeze(var_h(ss,:,:));    ppl_h_h=squeeze(ppl_h(ss,:,:));
% [M,I]=max(reshape(ppl_h_h,[],1));
% phi_h=reshape(var_h_h,[],1);
% %  y_h(ss)=  (phi_h(I));
%  y_h(ss)= circ_median (phi_h,1);
%
% end
y_h=phase_val_del_cr(ind_ses)';
y_m1=phase_val_del_cr(ind_ses_m1)';
y_m2=phase_val_del_cr(ind_ses_m2)';


rose(y_h,bins_num)
hold on
compass(aa*nanmean(exp(1j*y_h)),'k')

compass(aa*nanmean(exp(1j*y_m1)),'b')
compass(aa*nanmean(exp(1j*y_m2)),'g')

ax2=subplot(2,1,2);

var_wr=squeeze(plv_wr(ind_ses,1,:,:));
var_wr_m1=squeeze(plv_wr(ind_ses_m1,1,:,:));
var_wr_m2=squeeze(plv_wr(ind_ses_m2,1,:,:));


x_h=phase_val_del_wr(ind_ses)';
x_m1=phase_val_del_wr(ind_ses_m1)';
x_m2=phase_val_del_wr(ind_ses_m2)';

%     x_m1=(squeeze(circ_mean(circ_mean(var_wr_m1(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
%     x_m2=(squeeze(circ_mean(circ_mean(var_wr_m2(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
%
%     x_h=(squeeze(circ_mean(circ_mean(var_wr(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));

%     x_m1=(squeeze((circ_mean(var_wr_m1(:,ff,t1<t_h&t_h<t2),[],3))));
%     x_m2=(squeeze((circ_mean(var_wr_m2(:,ff,t1<t_h&t_h<t2),[],3))));
%
%     x_h=(squeeze((circ_mean(var_wr(:,ff,t1<t_h&t_h<t2),[],3))));
rose(x_h,bins_num)
hold on
compass(aa*nanmean(exp(1j*x_h)),'k')
compass(aa*nanmean(exp(1j*x_m1)),'b')
compass(aa*nanmean(exp(1j*x_m2)),'g')

%compass((exp(1j*x_h)))
axis([ax1 ax2],[bb1_ bb2_ bb1_ bb2_])
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],x_h,y_h],'fig3a',1)

%% von miss test
%
% y_h=[];y_h=(squeeze(circ_mean(circ_mean(var_cr(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
% y_h_m1=[];y_h_m1=(squeeze(circ_mean(circ_mean(var_cr_m1(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
% y_h_m2=[];y_h_m2=(squeeze(circ_mean(circ_mean(var_cr_m2(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));

% x_h=(squeeze(circ_mean(circ_mean(var_cr(ind_h,f>f1&f<f2,t_h<100),[],3),[],2)));

[pval_cr v_cr] = circ_vtest(y_h, circ_mean(y_h));
pval_cr
% [pval_cr z_cr] = circ_rtest(y_h)
kappa_cr_del_ss = circ_kappa(y_h)
phasediff_cr=circ_mean(y_h)*180/pi

kappa_cr_ss_m1 = circ_kappa(y_m1);
phasediff_cr_m1=circ_mean(y_m1)*180/pi;
kappa_cr_ss_m2 = circ_kappa(y_m2);
phasediff_cr_m2=circ_mean(y_m2)*180/pi;
%
% x_h=[];x_h=(squeeze(circ_mean(circ_mean(var_wr(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
% % x_h=(squeeze(circ_mean(circ_mean(var_wr(ind_h,f>f1&f<f2,t_h<100),[],3),[],2)));
% x_h_m1=[];x_h_m1=(squeeze(circ_mean(circ_mean(var_wr_m1(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));
% x_h_m2=[];x_h_m2=(squeeze(circ_mean(circ_mean(var_wr_m2(:,f>f1&f<f2,t1_delay<t_h&t_h<t2_delay),[],3),[],2)));

[pval_wr v_wr] = circ_vtest(x_h, circ_mean(x_h))
% [pval_wr z_wr] = circ_rtest(x_h)
kappa_wr_del_ss = circ_kappa(x_h)
phasediff_wr=circ_mean(x_h)*180/pi


kappa_wr_ss_m1 = circ_kappa(x_m1);
phasediff_wr_m1=circ_mean(x_m1)*180/pi;
kappa_wr_ss_m2 = circ_kappa(x_m2);
phasediff_wr_m2=circ_mean(x_m2)*180/pi;

diff_t=(mean(y_h)-mean(x_h))/sqrt((std(x_h)^2+std(y_h)^2)/50);
m_cr=mean(y_h);m_wr=mean(x_h);m_all=mean([x_h;y_h]);

%% Bootstrap test fo Kappa across sessions
y_base=phase_val_base_cr(ind_ses);
y_vis=phase_val_vis_cr(ind_ses);
y_del=phase_val_del_cr(ind_ses);

x_base=phase_val_base_wr(ind_ses);
x_vis=phase_val_vis_wr(ind_ses);
x_del=phase_val_del_wr(ind_ses);

N=length(y_del);

kappa_cr_base_ss=[];kappa_wr_base_ss=[];
kappa_cr_vis_ss=[];kappa_wr_vis_ss=[];
kappa_cr_del_ss=[];kappa_wr_del_ss=[];


for iteri=1:1001
    ind_h=randi(N,N,1);  
    
    kappa_cr_del_ss(iteri) = circ_kappa(y_del(ind_h));%-circ_kappa(x_0);
    kappa_cr_vis_ss(iteri) = circ_kappa(y_vis(ind_h));%-circ_kappa(x_0);
    kappa_cr_base_ss(iteri) = circ_kappa(y_base(ind_h));%-circ_kappa(x_0);
       
    kappa_wr_del_ss(iteri) = circ_kappa(x_del(ind_h));%-circ_kappa(x_0);
    kappa_wr_vis_ss(iteri) = circ_kappa(x_vis(ind_h));%-circ_kappa(x_0);
    kappa_wr_base_ss(iteri) = circ_kappa(x_base(ind_h));%-circ_kappa(x_0);
    
   
end

m_wr=nanmean(kappa_wr_del_ss);

m_cr=nanmean(kappa_cr_del_ss);

sum(kappa_cr_del_ss<m_wr)/length(kappa_cr_del_ss)

sum(kappa_wr_del_ss>m_cr)/length(kappa_cr_del_ss)

%% Fig. 3b Reliability of phase difference values across sessions
figure('name','Fig. 3b, Rezayat.E, et al')

ax=subplot(1,1,1);
hold on
errorbar([1 1.05 1.1],[nanmean(kappa_cr_base_ss) nanmean(kappa_cr_vis_ss) nanmean(kappa_cr_del_ss)],...
    [nanstd(kappa_cr_base_ss) nanstd(kappa_cr_vis_ss) nanstd(kappa_cr_del_ss)],'x','Color',Color_bhv(1,:))
errorbar([1 1.05 1.1],[nanmean(kappa_wr_base_ss) nanmean(kappa_wr_vis_ss) nanmean(kappa_wr_del_ss)],...
    [nanstd(kappa_wr_base_ss) nanstd(kappa_wr_vis_ss) nanstd(kappa_wr_del_ss)],'s','Color',Color_bhv(2,:))
ax.YTick=[0 1.4];
ax.XTickLabelRotation=30;
ax.XTick=[1 1.05 1.1];
ax.XTickLabel ={'Fixation','Visual','Delay'};
set(gca,'fontsize',14,'fontweight','bold')
ylabel('\kappa (a.u.)')
xlim([0.95 1.15])

export_file_to_excel([kappa_cr_base_ss' kappa_wr_base_ss' kappa_cr_vis_ss' kappa_wr_vis_ss' kappa_cr_del_ss' kappa_wr_del_ss'],'fig3b',12)

text(1.05,1.2,'Cr','color',Color_bhv(1,:),'fontsize',20)
text(1.05,1,'Wr','color',Color_bhv(2,:),'fontsize',20)
