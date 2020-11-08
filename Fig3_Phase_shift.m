function Fig3_Phase_shift(resp)
%% In the name of Allah
%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 3
% clear all
% close all
% clc

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
%% Load
gen_lfp_pairs(1).IT=resp(1).LFP_IT;
gen_lfp_pairs(1).FEF=resp(1).LFP_FEF;
gen_lfp_pairs(1).bhv{1}=resp(1).bhv;gen_lfp_pairs(1).condition{1}=resp(1).condition';
gen_lfp_pairs(1).preff_obj=resp(1).pref_obj_mu;
gen_lfp_pairs(1).npreff_obj=resp(1).npref_obj_mu;
gen_lfp_pairs(1).FEF_artifact=resp(1).FEF_artifact;
gen_lfp_pairs(1).IT_artifact=resp(1).IT_artifact;

%% Core calculation
disp('Wait! Calculating Phase shift  **********')
%%% shuffle within conditions
cal=1;
if cal
    ind_b=[50 110];                                                         % base line when down sample fs to 200 ( -300 -50)ms
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')
%     freqs2use=[4:50 51:3:100];
    freqs2use  =[4:1:70 71:5:130];                                         % frequency range
%    freqs2use  =logspace(log10(1),log10(100),60);
    st_ind=1; end_ind= 2100;                                               % time of task when fixatio is 500 ms
    
    %%% Definition of variabels
    plv_cr=[];plv_wr=[];                                                   % raw PLV for Cr& Wr
    plv_cr_z=[];plv_wr_z=[];                                               % baseline nomalized raw PLV for Cr& Wr
    plv_cr_n=[];plv_wr_n=[];                                               % shuffle corrected PLV for Cr& Wr
    plv_cr_n_z=[];plv_wr_n_z=[];                                           % shuffle corrected baseline nomalized PLV for Cr& Wr
    plv_shuff_cr_m=[];plv_shuff_wr_m=[];                                   % mean of shuffled  PLV for Cr& Wr
    plv_shuff_cr_std=[];plv_shuff_wr_std=[];                               % mean of shuffled baesline normalied PLV for Cr& Wr
    
    for ss=1:size(gen_lfp_pairs,2)
        
        
        tic
        %%% Determining  Condtions index
        lfp=gen_lfp_pairs(ss);
        ix_h_cr=[]; ix_h_wr=[];
        obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        for condi=1:6
            if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
            ix_h_cr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==1)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            ix_h_wr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==0)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            
        end
        
        
        %%% calculate morlet wavelt transform and spaontanous phase values
        %%% for each channel
        
        s_it_h = []; s_it_h = lfp.IT(:,st_ind:end_ind);
        s_fef_h = []; s_fef_h = lfp.FEF(:,st_ind:end_ind);
        p_it = []; p_fef = [];
        [p_it,p_fef,freqs2use_,t]= wavelet_core_DS_pertrial_phasedist(s_it_h,s_fef_h,freqs2use);
        phi_it = [];phi_it=(p_it);%phi_it=unwrap(phi_it);
        phi_fef =[];phi_fef=(p_fef);% phi_fef=unwrap(phi_fef);   
        
%   pow_it = []; pow_fef = []; phi_it = []; phi2_w = [];
%       [pow_it,pow_fef,phi_it,phi_fef,freqs2use,t]= wavelet_core_DS(s_it_h,s_fef_h);
%             

        N_trials=size(phi_it,1);
        
        %%% calculate PLV
        
        plv_cr_h=[];
        plv_cr_n_h=[];
        plv_shuff_cr_m_h=[];
        plv_shuff_cr_std_h=[];
        plv_wr_h=[];
        plv_wr_n_h=[];
        plv_shuff_wr_m_h=[];
        plv_shuff_wr_std_h=[];
        
        for condi=1:6
            N_trials_cr = [];N_trials_cr =  sum(logical(ix_h_cr(condi,:)));
            N_trials_wr = [];N_trials_wr =  sum(logical(ix_h_wr(condi,:)));
            
            N_sam=3;%min(N_trials_cr,N_trials_wr);
           % if(N_sam<3);N_sam=3;end
            
            phi_it_cr = [];phi_it_cr = squeeze(phi_it(logical(ix_h_cr(condi,:)),:,:));
            phi_fef_cr = [];phi_fef_cr = squeeze(phi_fef(logical(ix_h_cr(condi,:)),:,:));
            
            if N_trials_wr>0
                phi_it_wr = [];phi_it_wr = squeeze(phi_it(logical(ix_h_wr(condi,:)),:,:));
                phi_fef_wr = [];phi_fef_wr = squeeze(phi_fef(logical(ix_h_wr(condi,:)),:,:));
            end
            
            %%% Core PLV calculation over trials
            plv_h_cr=[]; plv_h_wr=[];            
            
            plv_h_cr=squeeze(angle(nanmean(exp(1j*((phi_fef_cr-phi_it_cr))),1)));
            
            if N_trials_wr>=N_sam
                plv_h_wr=squeeze(angle(nanmean(exp(1j*((phi_fef_wr-phi_it_wr))),1)));
             

            else
                 plv_h_wr=nan*ones(size(freqs2use,2),length(t));
                  plv_h_wr=nan*ones(size(freqs2use,2),length(t));
            end
            kappa_cr_h=[]; kappa_wr_h=[];
            for fi=1:size(phi_fef_cr,2)
                for ti=1:size(phi_fef_cr,3)
                    kappa_cr_h(fi,ti) = circ_kappa(angle(exp(1j*((phi_fef_cr(:,fi,ti)-phi_it_cr(:,fi,ti))))));
                    if N_trials_wr>=N_sam
                        kappa_wr_h(fi,ti) = circ_kappa(angle(exp(1j*((phi_fef_wr(:,fi,ti)-phi_it_wr(:,fi,ti))))));
                    else
                        kappa_wr_h(fi,ti)=nan;
                    end
                end
            end
            
            N_shuff=1001;
            plv_h_shuff_cr=[];plv_h_shuff_wr=[];
             
            %%% fill temporary variables and baseline normalization
            
            plv_cr_h(condi,:,:) = plv_h_cr;
            kappa_cr_h_(condi,:,:) = kappa_cr_h;

            
            plv_wr_h(condi,:,:) = plv_h_wr;
                kappa_wr_h_(condi,:,:) = kappa_wr_h;


        end
        
        %%% fill variables
        plv_cr(ss,:,:,:)=plv_cr_h;
                kappa_cr(ss,:,:,:)=kappa_cr_h_;

%        
        plv_wr(ss,:,:,:)=plv_wr_h;
        kappa_wr(ss,:,:,:)=kappa_wr_h_;

%         
        toc
        
    end
    
%     save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\phase_dist\phase_dist_full.mat'),'freqs2use',...
%         'plv_cr','plv_wr','kappa_cr','kappa_wr',...
%         'plv_shuff_cr_m','plv_shuff_cr_std','plv_shuff_wr_m','plv_shuff_wr_std')
    
end

%%
addpath(genpath('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\Codes\function\circstat-matlab-master'))
% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\phase_dist\phase_dist_full.mat'))
phase_cr=plv_cr;phase_wr=plv_wr;
% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\PPL_full_shuff_1001.mat'));
% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\PPL\PPL_shuff_1001.mat'));
load(strcat('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\PPL\PPL_shuff_1001_2.mat'));

% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\phase_dist\phase_dist_wavelet_over_trials.mat')

%% session selection
status = 2;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=1;%session_selection(4);
ind_ses_m2=setdiff(ind_ses,[1:51]);
ind_ses_m1=setdiff(ind_ses,[52:92]);
%%
ppl_cr=plv_cr;ppl_wr=plv_wr;
plv_cr=phase_cr;plv_wr=phase_wr;
NS=size(gen_lfp_pairs,2);
% Baseline Correction
baseline_corrected = 1;
ind_b = [41 110];
if baseline_corrected
    ppl_cr_h = ppl_cr;ppl_wr_h = ppl_wr;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
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

for ss=1:size(gen_lfp_pairs,2)
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

%% von miss test
try 
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
catch
end
%% Bootstrap test fo Kappa across sessions
try
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
catch
end
%% Jack knife test fo Kappa across sessions
if false
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


for iteri=1:length(y_base)
    ind_h=setdiff([1:length(y_base)],iteri);  
    
    kappa_cr_del_ss(iteri) = circ_kappa(y_del(ind_h));%-circ_kappa(x_0);
    kappa_cr_vis_ss(iteri) = circ_kappa(y_vis(ind_h));%-circ_kappa(x_0);
    kappa_cr_base_ss(iteri) = circ_kappa(y_base(ind_h));%-circ_kappa(x_0);
       
    kappa_wr_del_ss(iteri) = circ_kappa(x_del(ind_h));%-circ_kappa(x_0);
    kappa_wr_vis_ss(iteri) = circ_kappa(x_vis(ind_h));%-circ_kappa(x_0);
    kappa_wr_base_ss(iteri) = circ_kappa(x_base(ind_h));%-circ_kappa(x_0);
    
   
end

m_wr=nanmean(kappa_wr_del_ss);
prctile(kappa_cr_del_ss,2.5)
prctile(kappa_cr_del_ss,97.5)

prctile(kappa_wr_del_ss,2.5)
prctile(kappa_wr_del_ss,97.5)


m_cr=nanmean(kappa_cr_del_ss);
std_cr=nanstd(kappa_cr_del_ss);
m_cr-2*std_cr
m_cr+2*std_cr
sum(kappa_cr_del_ss<m_wr)/length(kappa_cr_del_ss)

sum(kappa_wr_del_ss>m_cr)/length(kappa_cr_del_ss)
end

%% Fig. 3b Reliability of phase difference values across sessions

figure('name','Fig. 3b, Rezayat.E, et al')
try 
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

text(1.05,1.2,'Cr','color',Color_bhv(1,:),'fontsize',20)
text(1.05,1,'Wr','color',Color_bhv(2,:),'fontsize',20)
catch 
end