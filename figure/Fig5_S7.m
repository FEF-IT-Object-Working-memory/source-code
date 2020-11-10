%% In the name of Allah
%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 2, Fig S3
clear all
close all
clc
%% Settings
ind_b=[1 500];
win_=300;step=10;
Percenti=33.33;
control=0;
Shuffle_highlow_PPL=0;
Shuffle_Preferance=1;
baseline_corrected_plv=1;
N_shuff=1001;
inf_cal=@(x,y,z) ROC_(x,y); method='roc';inf_pre=@(x,y) (x-y)./(y+x);

% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_wideband_full_test.mat')
% ind_ses=session_selection(1);

%% Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')


normalize=0;
mpi=1;

load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
%      load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))

load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_interval_bts_val',method,'_',num2str(Percenti),'.mat'))
bootstrap_in=(mi_high_in_bts-mi_low_in_bts)./(mi_high_in_bts+mi_low_in_bts);
bootstrap_out=(mi_high_out_bts-mi_low_out_bts)./(mi_high_out_bts+mi_low_out_bts);


mi_all_in_d= mi_all_in; mi_all_out_d= mi_all_out; mi_all_both_d= mi_all_both;
mi_high_in_d= mi_high_in; mi_high_out_d= mi_high_out; mi_high_both_d= mi_high_both;
mi_low_in_d= mi_low_in; mi_low_out_d= mi_low_out; mi_low_both_d= mi_low_both;

mi_all_in_d_n= (mi_all_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_all_out_d_n= (mi_all_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_all_both_d_n= (mi_all_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_high_in_d_n= (mi_high_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_high_out_d_n= (mi_high_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_high_both_d_n= (mi_high_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_low_in_d_n= (mi_low_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_low_out_d_n= (mi_low_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_low_both_d_n= (mi_low_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);


ind_pval=ones(253,1);
load([ cd ,'\Mat\PPL\ROC_interval_300_10ms_wideband_full_shuffled.mat'])


load(strcat([ cd ,'\Mat\PPL\ROC_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat']))
load(strcat([ cd ,'\Mat\PPL\Rate_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat']))

load([ cd ,'\Mat\PPL\plv_categorized_it_rate_wavelet_wideband_full.mat'])

%%Low number trials
in_rf=[];out_rf=[];
trail_num_in=[];trail_num_out=[];
for pp=1:size(val_it_all,1)
    val_in=[]; val_in=[ size(val_it_high{pp,1},1),size(val_it_low{pp,1},1),size(val_it_high{pp,3},1),size(val_it_low{pp,3},1)];
    val_out=[]; val_out=[  size(val_it_high{pp,4},1),size(val_it_low{pp,4},1),size(val_it_high{pp,6},1),size(val_it_low{pp,6},1); ];
    trail_num_in(pp)=min(min(val_in));
    trail_num_out(pp)=min(min(val_out));
    in_rf(pp)=(signrank(nanmean(val_it_all{pp,1}(:,1:500),2),nanmean(val_it_all{pp,1}(:,500:1000),2))<0.05);
    out_rf(pp)=(signrank(nanmean(val_it_all{pp,4}(:,1:500),2),nanmean(val_it_all{pp,4}(:,500:1000),2))<0.05);
end

ix_h_tri_in=(trail_num_in>1)';
ix_h_tri_out=(trail_num_out>1)';
% ix_h_tri_in=ones(253,1); ix_h_tri_out=ones(253,1);

%%Session selection
interval=1;
M1=0;
M2=0;
tresh=-1;
status = 1;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=session_selection(3);
%    ind_ses=setdiff(ind_ses,[74:92]);

ind_ses_m2=setdiff(ind_ses,[1:51]);
ind_ses_m1=setdiff(ind_ses,[52:92]);

if M1
    ind_ses=ind_ses_m1;
elseif M2
    ind_ses=ind_ses_m2;
    
end

%%% Remove nonselective LFP (same  to fig 2)
load multiunit_preferance_pval.mat
p=0.05;
ind_nonslective=find((p_rate>p)&(p_baseline>p));
ind_ses=setdiff(ind_ses,ind_nonslective);

t1=1100;t2=1700;
pp=0;ind_it_neurons_ses=[];ind_it_p=[];pval_dsi_it=[];
load([ cd,'\Data\index_it_neurons.mat'])

ind_it_neurons=ismember(ind_it_neurons_ses,ind_ses);

%%Neurons selection
% close all
clc
t1_delay=600;t2_delay=1200;
t1_selec=600;t2_selec=1200;

t_ind=find(t_h==751|t_h==1051);

if ~interval
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(nanmean(mi_all_both(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(nanmean(mi_all_in(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(nanmean(mi_all_out(:,t_ind),2)>=tresh)&ind_pval;
    
else
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(mi_all_both_d>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(mi_all_in_d>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(mi_all_out_d>=tresh)&ind_pval;
    
end

p_val=0.05;z_t=1.96;
p_val_in=ones(1,253);p_val_out=ones(1,253);p_val_both=ones(1,253);
for pp=1:253
    p_val_in(pp)=(sum(mi_all_in_d(pp)<squeeze(mi_in_shuff(pp,1,:)))/length(squeeze(mi_in_shuff(pp,1,:)))<p_val);
    p_val_out(pp)=(sum(mi_all_out_d(pp)<squeeze(mi_out_shuff(pp,1,:)))/length(squeeze(mi_out_shuff(pp,1,:)))<p_val);
    p_val_both(pp)=(sum(mi_all_both_d(pp)<squeeze(mi_both_shuff(pp,1,:)))/length(squeeze(mi_both_shuff(pp,1,:)))<p_val);
    
    % p_val_in(pp)=(mi_all_in_d(pp)>=prctile(squeeze(mi_in_shuff(pp,1,:)),90));
    % p_val_out(pp)=(mi_all_out_d(pp)>=prctile(squeeze(mi_out_shuff(pp,1,:)),90));
    % %
    % val_h=[];val_h= squeeze(mi_in_shuff(pp,1,:));ci_in=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_out_shuff(pp,1,:));ci_out=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_both_shuff(pp,1,:));ci_both=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % p_val_in(pp)=(mi_all_in_d(pp)>ci_in);
    % p_val_out(pp)=(mi_all_out_d(pp)>ci_out);
    % p_val_both(pp)=(mi_all_both_d(pp)>ci_both);
end

ix_resp=[];ix_resp=ones(length(ind_it_neurons),1);
th_h=0.002;p_val_d=0.05;
ih_s=(find(abs(nanmean(squeeze(rate_all_in(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
ih_s=(find(abs(nanmean(squeeze(rate_all_out(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
% ih_s=(find(abs(nanmean(squeeze(rate_all_both(:,1,t_h>=0&t_h<=1300)),2))<th_h));
% ix_resp(ih_s)=0;

ix_h_all=[];ix_h_all=ind_it_neurons'&ix_h_tri_in&ix_h_tri_out&ix_resp;
% ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&((p_val_in)'&(p_val_out)');
%ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&p_val_both';
ix_h_f=[];ix_h_f=p_val_both';

%     ix_h_f=[];ix_h_f=ix_h_i;
ix_h_in=ix_h_f&ix_h_all;
ix_h_out=ix_h_f&ix_h_all;
ix_h_both=ix_h_f&ix_h_all;
sum(ix_h_in)
%%%% Control for selecting one unit from each session
if control
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)'&ix_h_in);
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    
    %  for ui=1:length(ind_it_neurons_ses)
    %      ss=ind_it_neurons_ses(ui);
    %      if ix_h_in(ui)
    %          if ~ismember(ss,ses)
    %              ses=[ses,ss];
    %              idx_control(ui)=1;
    %
    %          end
    %      end
    %
    %  end
    ix_h_in=logical(idx_control)&ix_h_in;
    ix_h_out=logical(idx_control)&ix_h_in;
    ix_h_both=logical(idx_control)&ix_h_in;
    
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)');
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    ix_h_all=logical(idx_control)&ix_h_all;
end
%%Preview  ROC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preview scatter ROC

bin_num=20;
mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
mi_diff=mi_high_in-mi_low_in;
% mi_all_n_=nanmean(mi_all_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_low_n_=nanmean(mi_low_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_high_n_=nanmean(mi_high_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_diff_n_=nanmean(mi_diff(:,find(t_h>t1_delay&t_h<t2_delay)),2);
%
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
else
    
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
end

% signrank(mi_diff_n_)

% figure('name','ROC IN')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_in);
y_h_=[];y_h_=mi_high_n_(ix_h_in);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
% plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);

eval(strcat('selectivity_low_high_in= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))
%
% scatter(x_h,y_h,'k')
% text( 0.5,0.5,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,0.4,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')
%
% set(gca,'fontsize',14,'fontweight','bold')
% xlabel('AUC Low PPL')
% ylabel('AUC High PPL')
%  axis([0.2 1.1 0.2 1.1])
% % axis([-0.5 1.1 -0.5 1.1])

%  xbin=[-0.5:0.03:0.5];
% subplot(122)
% var_=[]; var_=-1*(y_h-x_h);
% % hist(var_,bin_num)
% bin_=0.05;
% nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
% [counts,centers]=hist(var_,nbin2);
% counts=counts/sum(counts);
% bar(centers,counts,1)
% line([0 0],ylim)
%
% m=median(var_);
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m m],ylim,'linewidth',2,'color','r')
% if(p<0.05)
%     text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
%
%     text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
% end
% set(gca,'fontsize',14,'fontweight','bold')
% ylabel('#')
% xlim([-0.4 0.4])

%%% OUT
if ~interval
    
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end
%
% figure('name','OUT')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_out);
y_h_=[];y_h_=mi_high_n_(ix_h_out);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);
eval(strcat('selectivity_low_high_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

%% Fig 5a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preveiw time course roc
bb_y_low=0.47;bb_y_high=0.85;

win=1;
figure ('name','Fig 5a')
ax_=subplot(1,1,1);
hold on
niceplot2(squeeze(mi_low_in(ix_h_in,:)),t_h'+150,win,0.7,0,0,':');
niceplot(squeeze(mi_high_in(ix_h_in,:)),t_h'+150,win,1,0.0,0);
bb1_=bb_y_low;bb2_=bb_y_high;
bb1=bb2_;
var_h_high=squeeze(mi_high_in(ix_h_in,:));
var_h_low=squeeze(mi_low_in(ix_h_in,:));
n_time=size(var_h_high,2);
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
% sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);

%    niceplot(squeeze(mi_all_in(ix_h_in,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
% preveiw time course OUT

export_file_to_excel([t_h nanmean(squeeze(mi_high_in(ix_h_in,:)))' nanmean(squeeze(mi_low_in(ix_h_in,:)))'],'fig5a',15)

%% Fig S7a

win=1;
figure ('name','Fig S7a')

ax_=subplot(1,1,1);
hold on
% plot(t_h,nanmean(mi_high_out(ix_h_out,:)),'g')
% plot(t_h,nanmean(mi_low_out(ix_h_out,:)),'b')
%   plot(t_h,nanmean(mi_all_out(ix_h_in,:)),'color',[0.3,0.3,0.3],'LineWidth',2)

% legend('High','Low')

niceplot2(squeeze(mi_low_out(ix_h_out,:)),t_h'+150,win,0,0,0.7,':');
niceplot(squeeze(mi_high_out(ix_h_out,:)),t_h'+150,win,0,0,1);
var_h_high=squeeze(mi_high_out(ix_h_out,:));
var_h_low=squeeze(mi_low_out(ix_h_out,:));
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;

%  niceplot(squeeze(mi_all_n(ix_h,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

% ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
export_file_to_excel([t_h nanmean(squeeze(mi_high_out(ix_h_out,:)))' nanmean(squeeze(mi_low_out(ix_h_out,:)))'],'figS7a',15)

%
bb_y_low=-0.1;bb_y_high=0.1;

%% Fig 5d
win=1;
figure ('name','Fig 5d')
mi_in_=inf_pre(mi_high_in(ix_h_both,:),mi_low_in(ix_h_both,:));
mi_out_=inf_pre(mi_high_out(ix_h_both,:),mi_low_out(ix_h_both,:));
mi_both_=inf_pre(mi_high_both(ix_h_both,:),mi_low_both(ix_h_both,:));
%
% [in_nan ~]=find(isnan(mi_in_));in_nan=unique(in_nan);
% [out_nan ~]=find(isnan(mi_out_));out_nan=unique(out_nan);
% mi_in_=mi_in_(setdiff([1:size(mi_in_,1)],[in_nan; out_nan] ),:);
% mi_out_=mi_out_(setdiff([1:size(mi_out_,1)],[in_nan; out_nan]),:);


hold on
% plot(t_h,nanmean(mi_in_),'r')
% plot(t_h,nanmean(mi_out_),'b')
%
%
% legend('IN','OUT')

niceplot(squeeze(mi_out_),t_h'+150,win,0,0,1);
niceplot(squeeze(mi_in_),t_h'+150,win,1,0,0);
export_file_to_excel([t_h nanmean(mi_in_)' nanmean(mi_out_)'],'fig5d',16)

% niceplot(squeeze(mi_both_),t_h',win,0.5,0.5,0.5);
bb1_=bb_y_low;bb2_=bb_y_high;
bb1=bb2_;

var_h_high=squeeze(mi_in_);
var_h_low=squeeze(mi_out_);
n_time=size(var_h_high,2);
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_high(:,tt),pp);end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)

line([-300 1500],[0 0],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];

ylabel('MI [High vs Low] (a.u)')
xlabel('Time from sample onset (ms)')

set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
%

if ~interval
    mi_out_n_=nanmean(mi_out_(:,t_ind),2);
    mi_in_n_=nanmean(mi_in_(:,t_ind),2);
else
    mi_in_n_=inf_pre(mi_high_in_d(ix_h_both,:),mi_low_in_d(ix_h_both,:));
    mi_out_n_=inf_pre(mi_high_out_d(ix_h_both,:),mi_low_out_d(ix_h_both,:));
    mi_both_n_=inf_pre(mi_high_both_d(ix_h_both,:),mi_low_both_d(ix_h_both,:));
end

m1_ind=((ind_it_neurons_ses(ix_h_both))<52);
m2_ind=((ind_it_neurons_ses(ix_h_both))>51);
%
%
% [p ,h,stats]=ranksum(mi_in_n_,mi_out_n_)
%
% nanmean(mi_in_n_)
% nanstd(mi_in_n_)/sqrt(length(mi_in_n_))
%
% nanmean(mi_out_n_)
% nanstd(mi_out_n_)/sqrt(length(mi_out_n_))
ind_neurons_resp=ind_it_neurons_ses(ix_h_both);
% save([ cd ,'PPL\control_selection.mat','mi_in_n_','mi_out_n_','ind_neurons_resp')
%% Fig 5e
figure ('name','Fig 5e')
subplot(1,1,1)
hold on
x_h_=[];x_h_=mi_out_n_;
y_h_=[];y_h_=mi_in_n_;
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
plot([-1:.1:1],[-1:0.1:1])
p=signrank(x_h,y_h)
n=length(y_h);
aaaa=sum(x_h<y_h)/n
%cohen_d1=(nanmean(y_h)-nanmean(x_h))/sqrt((((length(y_h)-1)*var(y_h))+((length(x_h)-1)*var(x_h)))/(length(y_h)+length(x_h)-2))
cohen_d=((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h)
var_in_bts_test=squeeze(bootstrap_in(ix_h_both,:,:));
var_out_bts_test=squeeze(bootstrap_out(ix_h_both,:,:));

var_in_bts_test=var_in_bts_test(~isnan(x_h_)&~isnan(y_h_),:);
var_out_bts_test=var_out_bts_test(~isnan(x_h_)&~isnan(y_h_),:);

for bts=1:size(var_in_bts_test,1)
    
    try
        % pval_right(bts)=signrank(var_in_bts_test(bts,:),var_out_bts_test(bts,:));
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        pval_right(bts)=(nanmean(var_in_bts_test(bts,:))>prctile(var_out_bts_test(bts,:),95));
        pval_left(bts)=(nanmean(var_out_bts_test(bts,:))>prctile(var_in_bts_test(bts,:),95));
        
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        
    catch
        pval_right(bts)=1;
        pval_left(bts)=1;
        
    end
end
% percent_of_increase=sum(pval_right<0.05)/length(pval_right)
% percent_of_decrease=sum(pval_left<0.05)/length(pval_right)
percent_of_increase=sum(pval_right==1)/length(pval_right)
percent_of_decrease=sum(pval_left==1)/length(pval_right)

eval(strcat('selectivity_mi_in_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')
export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5e',7)

text( 0.2,0.3,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.2,0.2,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')

set(gca,'fontsize',14,'fontweight','bold')
xlabel('MI OUT')
ylabel('MI IN')
axis([-0.4 0.6 -0.4 0.6])

figure ('name','Fig 5e inset')
subplot(1,1,1)
hold on
var_=[]; var_=-1*(y_h-x_h);
% hist(var_,bin_num)
bin_=0.05;
nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
[counts,centers]=hist(var_,nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)

xlim([-2 2.5])

m=median(var_);
line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
line([m m],ylim,'linewidth',2,'color','r')

m1=median(var_(m1_ind));
line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=median(var_(m2_ind));
line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
line([m2 m2],ylim,'linewidth',2,'color','g')

if(p<0.05)
    text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
    
    text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
end
set(gca,'fontsize',14,'fontweight','bold')
ylabel('% Pop.')
xlim([-0.5 0.5])

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
else
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
%  ix_h([9,178])=false;
%ix_h=ix_h_in;
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));

%% Fig 5f
figure('name','Fig 5f')
subplot(1,1,1)
title('IN')
hold on
% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')
export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5f',17)

[r1, m ,b]=regression(x_h,y_h,'one');
[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

% plotregression(x_h,y_h)
%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coecloficients{2,1};
%%%lm.Coefficients{1,1};

xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

length(x_h)
plot(xx,m*xx+b,'r')
plot(xx,m_m1*xx+b_m1,'b')
plot(xx,m_m2*xx+b_m2,'g')
xlabel(' IT selectivity (a.u)')
ylabel('Diff IT selectivity[High- Low](a.u)')
ylim([-0.3 0.6])
xlim([0.1 1])
% ylim([-0.4 0.4])
% xlim([0.2 1])
line([0.5 0.5],ylim)
line(xlim,[0 0])
set(gca,'fontsize',14,'fontweight','bold')
text( 0.5,-0.05,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.5,-0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.5,-0.15,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')

%% Fig S7d Correlation acrsos all It neurons OUt
figure('name','Fig S7d')

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);

x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
subplot(1,1,1)
% figure('name','OUT')
title('OUT')
hold on
% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')
export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'figS7d',17)

[r, m ,b]=regression(x_h,y_h,'one');

[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coefficients{2,1};
xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal');
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

length(x_h)

plot(xx,m*xx+b,'r')

plot(xx,m_m1*xx+b_m1,'b')
plot(xx,m_m2*xx+b_m2,'g')

xlabel(' IT selectivity (a.u)')
ylabel('Diff IT selectivity[High- Low](a.u)')
ylim([-0.3 0.6])
xlim([0.1 1])
line([0.5 0.5],ylim)
line(xlim,[0 0])
set(gca,'fontsize',14,'fontweight','bold')
text( 0.5,-0.05,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.5,-0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.5,-0.15,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')


%% IT power Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')


normalize=0;



load(strcat([ cd ,'\Mat\PPL\ROC_timeCoures_300_10ms_powit_full_interval_',method,'_',num2str(Percenti),'.mat']))
load(strcat([ cd ,'\Mat\PPL\ROC_shuff_powit_full_interval_',method,'_',num2str(Percenti),'.mat']))

mi_all_in_d= mi_all_in; mi_all_out_d= mi_all_out; mi_all_both_d= mi_all_both;
mi_high_in_d= mi_high_in; mi_high_out_d= mi_high_out; mi_high_both_d= mi_high_both;
mi_low_in_d= mi_low_in; mi_low_out_d= mi_low_out; mi_low_both_d= mi_low_both;

mi_all_in_d_n= (mi_all_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_all_out_d_n= (mi_all_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_all_both_d_n= (mi_all_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_high_in_d_n= (mi_high_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_high_out_d_n= (mi_high_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_high_both_d_n= (mi_high_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_low_in_d_n= (mi_low_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_low_out_d_n= (mi_low_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_low_both_d_n= (mi_low_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);


ind_pval=ones(253,1);
load([ cd ,'\Mat\PPL\ROC_interval_300_10ms_wideband_full_shuffled.mat'])


load(strcat([ cd ,'\Mat\PPL\ROC_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat']))
load(strcat([ cd ,'\Mat\PPL\Rate_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat']))

load([ cd ,'\Mat\PPL\plv_categorized_it_rate_wavelet_wideband_full.mat'])

%%Low number trials
in_rf=[];out_rf=[];
trail_num_in=[];trail_num_out=[];
for pp=1:size(val_it_all,1)
    val_in=[]; val_in=[ size(val_it_high{pp,1},1),size(val_it_low{pp,1},1),size(val_it_high{pp,3},1),size(val_it_low{pp,3},1)];
    val_out=[]; val_out=[  size(val_it_high{pp,4},1),size(val_it_low{pp,4},1),size(val_it_high{pp,6},1),size(val_it_low{pp,6},1); ];
    trail_num_in(pp)=min(min(val_in));
    trail_num_out(pp)=min(min(val_out));
    in_rf(pp)=(signrank(nanmean(val_it_all{pp,1}(:,1:500),2),nanmean(val_it_all{pp,1}(:,500:1000),2))<0.05);
    out_rf(pp)=(signrank(nanmean(val_it_all{pp,4}(:,1:500),2),nanmean(val_it_all{pp,4}(:,500:1000),2))<0.05);
end

ix_h_tri_in=(trail_num_in>1)';
ix_h_tri_out=(trail_num_out>1)';
% ix_h_tri_in=ones(253,1); ix_h_tri_out=ones(253,1);

%%Session selection
interval=1;
M1=0;
M2=0;
tresh=-1;
status = 1;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=session_selection(3);
%    ind_ses=setdiff(ind_ses,[74:92]);

ind_ses_m2=setdiff(ind_ses,[1:51]);
ind_ses_m1=setdiff(ind_ses,[52:92]);

if M1
    ind_ses=ind_ses_m1;
elseif M2
    ind_ses=ind_ses_m2;
    
end

%%% Remove nonselective LFP (same  to fig 2)
load multiunit_preferance_pval.mat
p=0.05;
ind_nonslective=find((p_rate>p)&(p_baseline>p));
ind_ses=setdiff(ind_ses,ind_nonslective);

%%Select based on sessions
t1=1100;t2=1700;

ind_it_neurons=ismember(ind_it_neurons_ses,ind_ses);

%%Neurons selection
% close all
clc
t1_delay=600;t2_delay=1200;
t1_selec=600;t2_selec=1200;

t_ind=find(t_h==751|t_h==1051);

if ~interval
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(nanmean(mi_all_both(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(nanmean(mi_all_in(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(nanmean(mi_all_out(:,t_ind),2)>=tresh)&ind_pval;
    
else
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(mi_all_both_d>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(mi_all_in_d>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(mi_all_out_d>=tresh)&ind_pval;
    
end

p_val=0.05;z_t=1.96;
p_val_in=ones(1,253);p_val_out=ones(1,253);p_val_both=ones(1,253);
for pp=1:253
    p_val_in(pp)=(sum(mi_all_in_d(pp)<squeeze(mi_in_shuff(pp,1,:)))/length(squeeze(mi_in_shuff(pp,1,:)))<p_val);
    p_val_out(pp)=(sum(mi_all_out_d(pp)<squeeze(mi_out_shuff(pp,1,:)))/length(squeeze(mi_out_shuff(pp,1,:)))<p_val);
    p_val_both(pp)=(sum(mi_all_both_d(pp)<squeeze(mi_both_shuff(pp,1,:)))/length(squeeze(mi_both_shuff(pp,1,:)))<p_val);
    
    % p_val_in(pp)=(mi_all_in_d(pp)>=prctile(squeeze(mi_in_shuff(pp,1,:)),90));
    % p_val_out(pp)=(mi_all_out_d(pp)>=prctile(squeeze(mi_out_shuff(pp,1,:)),90));
    % %
    % val_h=[];val_h= squeeze(mi_in_shuff(pp,1,:));ci_in=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_out_shuff(pp,1,:));ci_out=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_both_shuff(pp,1,:));ci_both=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % p_val_in(pp)=(mi_all_in_d(pp)>ci_in);
    % p_val_out(pp)=(mi_all_out_d(pp)>ci_out);
    % p_val_both(pp)=(mi_all_both_d(pp)>ci_both);
end

ix_resp=[];ix_resp=ones(length(ind_it_neurons),1);
th_h=0.002;p_val_d=0.05;
ih_s=(find(abs(nanmean(squeeze(rate_all_in(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
ih_s=(find(abs(nanmean(squeeze(rate_all_out(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
% ih_s=(find(abs(nanmean(squeeze(rate_all_both(:,1,t_h>=0&t_h<=1300)),2))<th_h));
% ix_resp(ih_s)=0;

ix_h_all=[];ix_h_all=ind_it_neurons'&ix_h_tri_in&ix_h_tri_out&ix_resp;
% ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&((p_val_in)'&(p_val_out)');
%ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&p_val_both';
ix_h_f=[];ix_h_f=p_val_both';

%     ix_h_f=[];ix_h_f=ix_h_i;
ix_h_in=ix_h_f&ix_h_all;
ix_h_out=ix_h_f&ix_h_all;
ix_h_both=ix_h_f&ix_h_all;
sum(ix_h_in)
%%%% Control for selecting one unit from each session
if control
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)'&ix_h_in);
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    
    %  for ui=1:length(ind_it_neurons_ses)
    %      ss=ind_it_neurons_ses(ui);
    %      if ix_h_in(ui)
    %          if ~ismember(ss,ses)
    %              ses=[ses,ss];
    %              idx_control(ui)=1;
    %
    %          end
    %      end
    %
    %  end
    ix_h_in=logical(idx_control)&ix_h_in;
    ix_h_out=logical(idx_control)&ix_h_in;
    ix_h_both=logical(idx_control)&ix_h_in;
    
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)');
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    ix_h_all=logical(idx_control)&ix_h_all;
end


bin_num=20;
mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
mi_diff=mi_high_in-mi_low_in;
% mi_all_n_=nanmean(mi_all_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_low_n_=nanmean(mi_low_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_high_n_=nanmean(mi_high_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_diff_n_=nanmean(mi_diff(:,find(t_h>t1_delay&t_h<t2_delay)),2);
%
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
else
    
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
end

% signrank(mi_diff_n_)
%
% figure('name','ROC IN')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_in);
y_h_=[];y_h_=mi_high_n_(ix_h_in);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
% plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);

eval(strcat('IT_Power_selectivity_low_high_in= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))


%%% OUT
if ~interval
    
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end

% figure('name','OUT')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_out);
y_h_=[];y_h_=mi_high_n_(ix_h_out);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
% plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);
eval(strcat('selectivity_low_high_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))


%% Fig 5b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preveiw time course roc
bb_y_low=0.47;bb_y_high=0.85;

win=1;
figure ('name','Fig 5b')
ax_=subplot(1,1,1);
hold on
niceplot2(squeeze(mi_low_in(ix_h_in,:)),t_h'+150,win,0.7,0,0,':');
niceplot(squeeze(mi_high_in(ix_h_in,:)),t_h'+150,win,1,0.0,0);
bb1_=bb_y_low;bb2_=bb_y_high;
bb1=bb2_;
var_h_high=squeeze(mi_high_in(ix_h_in,:));
var_h_low=squeeze(mi_low_in(ix_h_in,:));
n_time=size(var_h_high,2);
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
% sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);
export_file_to_excel([t_h nanmean(squeeze(mi_high_in(ix_h_in,:)))' nanmean(squeeze(mi_low_in(ix_h_in,:)))'],'fig5b',15)

%    niceplot(squeeze(mi_all_in(ix_h_in,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
% preveiw time course OUT

%% Fig S7b

win=1;
figure('name','Fig S7b')

ax_=subplot(1,1,1);
hold on
% plot(t_h,nanmean(mi_high_out(ix_h_out,:)),'g')
% plot(t_h,nanmean(mi_low_out(ix_h_out,:)),'b')
%   plot(t_h,nanmean(mi_all_out(ix_h_in,:)),'color',[0.3,0.3,0.3],'LineWidth',2)

% legend('High','Low')

niceplot2(squeeze(mi_low_out(ix_h_out,:)),t_h'+150,win,0,0,0.7,':');
niceplot(squeeze(mi_high_out(ix_h_out,:)),t_h'+150,win,0,0,1);
var_h_high=squeeze(mi_high_out(ix_h_out,:));
var_h_low=squeeze(mi_low_out(ix_h_out,:));
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
export_file_to_excel([t_h nanmean(squeeze(mi_high_out(ix_h_out,:)))' nanmean(squeeze(mi_low_out(ix_h_out,:)))'],'figS7b',15)

%  niceplot(squeeze(mi_all_n(ix_h,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

% ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MI time coures
bb_y_low=-0.1;bb_y_high=0.1;


win=1;
mi_in_=inf_pre(mi_high_in(ix_h_both,:),mi_low_in(ix_h_both,:));
mi_out_=inf_pre(mi_high_out(ix_h_both,:),mi_low_out(ix_h_both,:));
mi_both_=inf_pre(mi_high_both(ix_h_both,:),mi_low_both(ix_h_both,:));
%
% [in_nan ~]=find(isnan(mi_in_));in_nan=unique(in_nan);
% [out_nan ~]=find(isnan(mi_out_));out_nan=unique(out_nan);
% mi_in_=mi_in_(setdiff([1:size(mi_in_,1)],[in_nan; out_nan] ),:);
% mi_out_=mi_out_(setdiff([1:size(mi_out_,1)],[in_nan; out_nan]),:);



var_h_high=squeeze(mi_in_);
var_h_low=squeeze(mi_out_);
n_time=size(var_h_high,2);
% p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_high(:,tt),pp);end
% sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
% sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);
% ylim([bb_y_low bb_y_high])
% xlim([-150 1550])
% line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)

% line([-300 1500],[0 0],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];

% ylabel('MI [High vs Low] (a.u)')
% xlabel('Time from sample onset (ms)')
%
% set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
%

if ~interval
    mi_out_n_=nanmean(mi_out_(:,t_ind),2);
    mi_in_n_=nanmean(mi_in_(:,t_ind),2);
else
    mi_in_n_=inf_pre(mi_high_in_d(ix_h_both,:),mi_low_in_d(ix_h_both,:));
    mi_out_n_=inf_pre(mi_high_out_d(ix_h_both,:),mi_low_out_d(ix_h_both,:));
    mi_both_n_=inf_pre(mi_high_both_d(ix_h_both,:),mi_low_both_d(ix_h_both,:));
end

m1_ind=((ind_it_neurons_ses(ix_h_both))<52);
m2_ind=((ind_it_neurons_ses(ix_h_both))>51);
%
%
% [p ,h,stats]=ranksum(mi_in_n_,mi_out_n_)
%
% nanmean(mi_in_n_)
% nanstd(mi_in_n_)/sqrt(length(mi_in_n_))
%
% nanmean(mi_out_n_)
% nanstd(mi_out_n_)/sqrt(length(mi_out_n_))
ind_neurons_resp=ind_it_neurons_ses(ix_h_both);
% save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\control_selection.mat','mi_in_n_','mi_out_n_','ind_neurons_resp')
% figure
% subplot(121)
% hold on
x_h_=[];x_h_=mi_out_n_;
y_h_=[];y_h_=mi_in_n_;
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
% plot([-1:.1:1],[-1:0.1:1])
p=signrank(x_h,y_h)
n=length(y_h);
aaaa=sum(x_h<y_h)/n
%cohen_d1=(nanmean(y_h)-nanmean(x_h))/sqrt((((length(y_h)-1)*var(y_h))+((length(x_h)-1)*var(x_h)))/(length(y_h)+length(x_h)-2))
cohen_d=((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h)
var_in_bts_test=squeeze(bootstrap_in(ix_h_both,:,:));
var_out_bts_test=squeeze(bootstrap_out(ix_h_both,:,:));

var_in_bts_test=var_in_bts_test(~isnan(x_h_)&~isnan(y_h_),:);
var_out_bts_test=var_out_bts_test(~isnan(x_h_)&~isnan(y_h_),:);

for bts=1:size(var_in_bts_test,1)
    
    try
        % pval_right(bts)=signrank(var_in_bts_test(bts,:),var_out_bts_test(bts,:));
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        pval_right(bts)=(nanmean(var_in_bts_test(bts,:))>prctile(var_out_bts_test(bts,:),95));
        pval_left(bts)=(nanmean(var_out_bts_test(bts,:))>prctile(var_in_bts_test(bts,:),95));
        
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        
    catch
        pval_right(bts)=1;
        pval_left(bts)=1;
        
    end
end
% percent_of_increase=sum(pval_right<0.05)/length(pval_right)
% percent_of_decrease=sum(pval_left<0.05)/length(pval_right)
percent_of_increase=sum(pval_right==1)/length(pval_right)
percent_of_decrease=sum(pval_left==1)/length(pval_right)

eval(strcat('selectivity_mi_in_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

% scatter(x_h,y_h)
% scatter(x_h(m1_ind),y_h(m1_ind),'b*')
% scatter(x_h(m2_ind),y_h(m2_ind),'gs')
export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5e',7)

text( 0.2,0.3,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.2,0.2,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')

% set(gca,'fontsize',14,'fontweight','bold')
% xlabel('MI OUT')
% ylabel('MI IN')
%  axis([-0.4 0.6 -0.4 0.6])
% subplot(122)
% var_=[]; var_=-1*(y_h-x_h);
% % hist(var_,bin_num)
% bin_=0.05;
% nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
% [counts,centers]=hist(var_,nbin2);
% counts=counts/sum(counts);
% bar(centers,counts,1)
% line([0 0],ylim)
%
% xlim([-2 2.5])
%
% m=median(var_);
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m m],ylim,'linewidth',2,'color','r')
%
% m1=median(var_(m1_ind));
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m1 m1],ylim,'linewidth',2,'color','b')
%
% m2=median(var_(m2_ind));
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m2 m2],ylim,'linewidth',2,'color','g')
%
% if(p<0.05)
%     text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
%
%     text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
% end
% set(gca,'fontsize',14,'fontweight','bold')
% ylabel('% Pop.')
% xlim([-0.5 0.5])

%%Fig behrad

%%Correlation acrsos all IT neurons IN

%  ix_h_all=[];ix_h_all=ind_it_neurons'&ix_h_tri_in&ix_h_tri_out&ind_pval;

% [m ~]=find(isnan(mi_high_out(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_low_out(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_high_in(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_low_in(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% %

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
else
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
%  ix_h([9,178])=false;
%ix_h=ix_h_in;
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));

% figure('name','IN')
% subplot(121)
% title('IN')
% hold on
% % scatter(x_h,y_h)
% scatter(x_h(m1_ind),y_h(m1_ind),'b*')
% scatter(x_h(m2_ind),y_h(m2_ind),'gs')
% export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5f',7)

[r1, m ,b]=regression(x_h,y_h,'one');
[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

% plotregression(x_h,y_h)
%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coecloficients{2,1};
%%%lm.Coefficients{1,1};

xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

% length(x_h)
% plot(xx,m*xx+b,'r')
% plot(xx,m_m1*xx+b_m1,'b')
% plot(xx,m_m2*xx+b_m2,'g')
% xlabel(' IT selectivity (a.u)')
% ylabel('Diff IT selectivity[High- Low](a.u)')
% ylim([-0.3 0.6])
%  xlim([0.1 1])
% % ylim([-0.4 0.4])
% % xlim([0.2 1])
% line([0.5 0.5],ylim)
% line(xlim,[0 0])
% set(gca,'fontsize',14,'fontweight','bold')
% text( 0.5,-0.05,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.15,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')


%%% Correlation acrsos all It neurons OUt

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);

x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
% subplot(122)
% % figure('name','OUT')
% title('OUT')
% hold on
% % scatter(x_h,y_h)
% scatter(x_h(m1_ind),y_h(m1_ind),'b*')
% scatter(x_h(m2_ind),y_h(m2_ind),'gs')
% export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'figS7d',7)

[r, m ,b]=regression(x_h,y_h,'one');

[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coefficients{2,1};
xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal');
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

length(x_h)

% plot(xx,m*xx+b,'r')
%
% plot(xx,m_m1*xx+b_m1,'b')
% plot(xx,m_m2*xx+b_m2,'g')
%
% xlabel(' IT selectivity (a.u)')
% ylabel('Diff IT selectivity[High- Low](a.u)')
% ylim([-0.3 0.6])
%  xlim([0.1 1])
% line([0.5 0.5],ylim)
% line(xlim,[0 0])
% set(gca,'fontsize',14,'fontweight','bold')
% text( 0.5,-0.05,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.15,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')


%% FEF Power Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')


normalize=0;
mpi=3;


load(strcat([ cd ,'\Mat\PPL\ROC_timeCoures_300_10ms_powfef_full_interval_',method,'_',num2str(Percenti),'.mat']))
load(strcat([ cd ,'\Mat\PPL\ROC_shuff_powfef_full_interval_',method,'_',num2str(Percenti),'.mat']))

mi_all_in_d= mi_all_in; mi_all_out_d= mi_all_out; mi_all_both_d= mi_all_both;
mi_high_in_d= mi_high_in; mi_high_out_d= mi_high_out; mi_high_both_d= mi_high_both;
mi_low_in_d= mi_low_in; mi_low_out_d= mi_low_out; mi_low_both_d= mi_low_both;

mi_all_in_d_n= (mi_all_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_all_out_d_n= (mi_all_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_all_both_d_n= (mi_all_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_high_in_d_n= (mi_high_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_high_out_d_n= (mi_high_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_high_both_d_n= (mi_high_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);
mi_low_in_d_n= (mi_low_in-nanmean(squeeze(mi_in_shuff),2))./nanstd(squeeze(mi_in_shuff),[],2);
mi_low_out_d_n= (mi_low_out-nanmean(squeeze(mi_out_shuff),2))./nanstd(squeeze(mi_out_shuff),[],2);
mi_low_both_d_n= (mi_low_both-nanmean(squeeze(mi_both_shuff),2))./nanstd(squeeze(mi_both_shuff),[],2);


ind_pval=ones(253,1);
load([ cd ,'\Mat\PPL\ROC_interval_300_10ms_wideband_full_shuffled.mat'])


load(strcat([ cd ,'\Mat\PPL\ROC_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat']))
load(strcat([ cd ,'\Mat\PPL\Rate_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat']))
load([ cd ,'\Mat\PPL\plv_categorized_it_rate_wavelet_wideband_full.mat'])

%%Low number trials
in_rf=[];out_rf=[];
trail_num_in=[];trail_num_out=[];
for pp=1:size(val_it_all,1)
    val_in=[]; val_in=[ size(val_it_high{pp,1},1),size(val_it_low{pp,1},1),size(val_it_high{pp,3},1),size(val_it_low{pp,3},1)];
    val_out=[]; val_out=[  size(val_it_high{pp,4},1),size(val_it_low{pp,4},1),size(val_it_high{pp,6},1),size(val_it_low{pp,6},1); ];
    trail_num_in(pp)=min(min(val_in));
    trail_num_out(pp)=min(min(val_out));
    in_rf(pp)=(signrank(nanmean(val_it_all{pp,1}(:,1:500),2),nanmean(val_it_all{pp,1}(:,500:1000),2))<0.05);
    out_rf(pp)=(signrank(nanmean(val_it_all{pp,4}(:,1:500),2),nanmean(val_it_all{pp,4}(:,500:1000),2))<0.05);
end

ix_h_tri_in=(trail_num_in>1)';
ix_h_tri_out=(trail_num_out>1)';
% ix_h_tri_in=ones(253,1); ix_h_tri_out=ones(253,1);

% Session selection
interval=1;
M1=0;
M2=0;
tresh=-1;
status = 1;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=session_selection(3);
%    ind_ses=setdiff(ind_ses,[74:92]);

ind_ses_m2=setdiff(ind_ses,[1:51]);
ind_ses_m1=setdiff(ind_ses,[52:92]);

if M1
    ind_ses=ind_ses_m1;
elseif M2
    ind_ses=ind_ses_m2;
    
end

%%% Remove nonselective LFP (same  to fig 2)
load multiunit_preferance_pval.mat
p=0.05;
ind_nonslective=find((p_rate>p)&(p_baseline>p));
ind_ses=setdiff(ind_ses,ind_nonslective);

%%Select based on sessions
t1=1100;t2=1700;
ind_it_neurons=ismember(ind_it_neurons_ses,ind_ses);

%%Neurons selection
% close all
clc
t1_delay=600;t2_delay=1200;
t1_selec=600;t2_selec=1200;

t_ind=find(t_h==751|t_h==1051);

if ~interval
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(nanmean(mi_all_both(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(nanmean(mi_all_in(:,t_ind),2)>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(nanmean(mi_all_out(:,t_ind),2)>=tresh)&ind_pval;
    
else
    
    ix_h_b=[];ix_h_b=ones(length(ind_it_neurons),1);ix_h_b=(mi_all_both_d>=tresh)&ind_pval;
    ix_h_i=[];ix_h_i=ones(length(ind_it_neurons),1);ix_h_i=(mi_all_in_d>=tresh)&ind_pval;
    ix_h_o=[];ix_h_o=ones(length(ind_it_neurons),1);ix_h_o=(mi_all_out_d>=tresh)&ind_pval;
    
end

p_val=0.05;z_t=1.96;
p_val_in=ones(1,253);p_val_out=ones(1,253);p_val_both=ones(1,253);
for pp=1:253
    p_val_in(pp)=(sum(mi_all_in_d(pp)<squeeze(mi_in_shuff(pp,1,:)))/length(squeeze(mi_in_shuff(pp,1,:)))<p_val);
    p_val_out(pp)=(sum(mi_all_out_d(pp)<squeeze(mi_out_shuff(pp,1,:)))/length(squeeze(mi_out_shuff(pp,1,:)))<p_val);
    p_val_both(pp)=(sum(mi_all_both_d(pp)<squeeze(mi_both_shuff(pp,1,:)))/length(squeeze(mi_both_shuff(pp,1,:)))<p_val);
    
    % p_val_in(pp)=(mi_all_in_d(pp)>=prctile(squeeze(mi_in_shuff(pp,1,:)),90));
    % p_val_out(pp)=(mi_all_out_d(pp)>=prctile(squeeze(mi_out_shuff(pp,1,:)),90));
    % %
    % val_h=[];val_h= squeeze(mi_in_shuff(pp,1,:));ci_in=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_out_shuff(pp,1,:));ci_out=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % val_h=[];val_h= squeeze(mi_both_shuff(pp,1,:));ci_both=nanmean(val_h)+z_t*nanstd(val_h);%/sqrt(length(val_h));
    % p_val_in(pp)=(mi_all_in_d(pp)>ci_in);
    % p_val_out(pp)=(mi_all_out_d(pp)>ci_out);
    % p_val_both(pp)=(mi_all_both_d(pp)>ci_both);
end

ix_resp=[];ix_resp=ones(length(ind_it_neurons),1);
th_h=0.002;p_val_d=0.05;
ih_s=(find(abs(nanmean(squeeze(rate_all_in(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
ih_s=(find(abs(nanmean(squeeze(rate_all_out(:,1,t_h>0&t_h<1300)),2))<th_h));
ix_resp(ih_s)=0;
% ih_s=(find(abs(nanmean(squeeze(rate_all_both(:,1,t_h>=0&t_h<=1300)),2))<th_h));
% ix_resp(ih_s)=0;

ix_h_all=[];ix_h_all=ind_it_neurons'&ix_h_tri_in&ix_h_tri_out&ix_resp;
% ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&((p_val_in)'&(p_val_out)');
%ix_h_f=[];ix_h_f=ix_h_i&ix_h_o&p_val_both';
ix_h_f=[];ix_h_f=p_val_both';

%     ix_h_f=[];ix_h_f=ix_h_i;
ix_h_in=ix_h_f&ix_h_all;
ix_h_out=ix_h_f&ix_h_all;
ix_h_both=ix_h_f&ix_h_all;
sum(ix_h_in)
%%%% Control for selecting one unit from each session
if control
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)'&ix_h_in);
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    
    %  for ui=1:length(ind_it_neurons_ses)
    %      ss=ind_it_neurons_ses(ui);
    %      if ix_h_in(ui)
    %          if ~ismember(ss,ses)
    %              ses=[ses,ss];
    %              idx_control(ui)=1;
    %
    %          end
    %      end
    %
    %  end
    ix_h_in=logical(idx_control)&ix_h_in;
    ix_h_out=logical(idx_control)&ix_h_in;
    ix_h_both=logical(idx_control)&ix_h_in;
    
    idx_control=zeros(length(ind_it_neurons_ses),1);
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_it_neurons_ses,ss)');
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    ix_h_all=logical(idx_control)&ix_h_all;
end
%%Preview  ROC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preview scatter ROC

bin_num=20;
mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
mi_diff=mi_high_in-mi_low_in;
% mi_all_n_=nanmean(mi_all_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_low_n_=nanmean(mi_low_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_high_n_=nanmean(mi_high_in(:,find(t_h>t1_delay&t_h<t2_delay)),2);
% mi_diff_n_=nanmean(mi_diff(:,find(t_h>t1_delay&t_h<t2_delay)),2);
%
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
else
    
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
    mi_diff_n_=nanmean(mi_diff(:,t_ind),2);
end

% signrank(mi_diff_n_)

% figure('name','ROC IN')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_in);
y_h_=[];y_h_=mi_high_n_(ix_h_in);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
% plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);

eval(strcat('selectivity_low_high_in= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

%
% %%% OUT
if ~interval
    
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end

% figure('name','OUT')
% subplot(121)
% hold on
x_h_=[];x_h_=mi_low_n_(ix_h_out);
y_h_=[];y_h_=mi_high_n_(ix_h_out);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
% plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);
eval(strcat('selectivity_low_high_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

%% Fig5c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preveiw time course roc
bb_y_low=0.47;bb_y_high=0.85;

win=1;
figure ('name','Fig5c')
ax_=subplot(1,1,1);
hold on
niceplot2(squeeze(mi_low_in(ix_h_in,:)),t_h'+150,win,0.7,0,0,':');
niceplot(squeeze(mi_high_in(ix_h_in,:)),t_h'+150,win,1,0.0,0);
bb1_=bb_y_low;bb2_=bb_y_high;
bb1=bb2_;
var_h_high=squeeze(mi_high_in(ix_h_in,:));
var_h_low=squeeze(mi_low_in(ix_h_in,:));
n_time=size(var_h_high,2);
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
% sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);
export_file_to_excel([t_h nanmean(squeeze(mi_high_in(ix_h_in,:)))' nanmean(squeeze(mi_low_in(ix_h_in,:)))'],'fig5c',15)

%    niceplot(squeeze(mi_all_in(ix_h_in,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
% preveiw time course OUT

%% FigS7c

win=1;
figure ('name','FigS7c')

ax_=subplot(1,1,1);
hold on
% plot(t_h,nanmean(mi_high_out(ix_h_out,:)),'g')
% plot(t_h,nanmean(mi_low_out(ix_h_out,:)),'b')
%   plot(t_h,nanmean(mi_all_out(ix_h_in,:)),'color',[0.3,0.3,0.3],'LineWidth',2)

% legend('High','Low')

niceplot2(squeeze(mi_low_out(ix_h_out,:)),t_h'+150,win,0,0,0.7,':');
niceplot(squeeze(mi_high_out(ix_h_out,:)),t_h'+150,win,0,0,1);
var_h_high=squeeze(mi_high_out(ix_h_out,:));
var_h_low=squeeze(mi_low_out(ix_h_out,:));
p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_low(:,tt),var_h_high(:,tt));end
sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
export_file_to_excel([t_h nanmean(squeeze(mi_high_out(ix_h_out,:)))' nanmean(squeeze(mi_low_out(ix_h_out,:)))'],'figS7c',15)

%  niceplot(squeeze(mi_all_n(ix_h,:)),t_h',win,1,0,0);
% text(900 ,nanmean(nanmean(mi_low_n(ix_h,:))),strcat('N=',num2str(sum(ix_h))),'color','r')
line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
line(xlim,[0.5 0.5],'color','k','LineWidth',2)

% ax_.Xtick=[0 300 1300];
ylabel('IT selectivity (a.u)')
xlabel('Time from sample onset (ms)')
ylim([bb_y_low bb_y_high])
xlim([-150 1550])
set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MI time coures
bb_y_low=-0.1;bb_y_high=0.1;


win=1;
% figure ('name','IN')
mi_in_=inf_pre(mi_high_in(ix_h_both,:),mi_low_in(ix_h_both,:));
mi_out_=inf_pre(mi_high_out(ix_h_both,:),mi_low_out(ix_h_both,:));
mi_both_=inf_pre(mi_high_both(ix_h_both,:),mi_low_both(ix_h_both,:));
%
% [in_nan ~]=find(isnan(mi_in_));in_nan=unique(in_nan);
% [out_nan ~]=find(isnan(mi_out_));out_nan=unique(out_nan);
% mi_in_=mi_in_(setdiff([1:size(mi_in_,1)],[in_nan; out_nan] ),:);
% mi_out_=mi_out_(setdiff([1:size(mi_out_,1)],[in_nan; out_nan]),:);


% plot(t_h,nanmean(mi_in_),'r')
% plot(t_h,nanmean(mi_out_),'b')
%
%
% legend('IN','OUT')
%
% niceplot(squeeze(mi_out_),t_h'+150,win,0,0,1);
% niceplot(squeeze(mi_in_),t_h'+150,win,1,0,0);
% % niceplot(squeeze(mi_both_),t_h',win,0.5,0.5,0.5);
% bb1_=bb_y_low;bb2_=bb_y_high;
% bb1=bb2_;

var_h_high=squeeze(mi_in_);
var_h_low=squeeze(mi_out_);
n_time=size(var_h_high,2);
% p=[];for tt=1:n_time;[p(tt) ~]=signrank(var_h_high(:,tt),pp);end
% sig_05=[];sig_05=nan*ones(1,n_time);sig_05(p<0.05)=1*bb1;
% sig_01=nan*ones(1,n_time);sig_01(p<0.01)=1*(bb1-(bb2_-bb1_)/10*1);
% % sig_001=nan*ones(1,n_time);sig_001(p<0.001)=1*(bb1-(bb2_-bb1_)/10*2);
% ylim([bb_y_low bb_y_high])
% xlim([-150 1550])
% line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
% line([600 1200],[bb_y_low bb_y_low],'color','k','LineWidth',4)
%
% line([-300 1500],[0 0],'color','k','LineWidth',2)

%  ax_.Xtick=[0 300 1300];

% ylabel('MI [High vs Low] (a.u)')
% xlabel('Time from sample onset (ms)')
%
% set(gca,'fontsize',14,'fontweight','bold')
% plot(t_h,sig_05,'*k')
% plot(t_h,sig_01,'+k')
% plot(t_h,sig_001,'sk')
%

if ~interval
    mi_out_n_=nanmean(mi_out_(:,t_ind),2);
    mi_in_n_=nanmean(mi_in_(:,t_ind),2);
else
    mi_in_n_=inf_pre(mi_high_in_d(ix_h_both,:),mi_low_in_d(ix_h_both,:));
    mi_out_n_=inf_pre(mi_high_out_d(ix_h_both,:),mi_low_out_d(ix_h_both,:));
    mi_both_n_=inf_pre(mi_high_both_d(ix_h_both,:),mi_low_both_d(ix_h_both,:));
end

m1_ind=((ind_it_neurons_ses(ix_h_both))<52);
m2_ind=((ind_it_neurons_ses(ix_h_both))>51);
%
%
% [p ,h,stats]=ranksum(mi_in_n_,mi_out_n_)
%
% nanmean(mi_in_n_)
% nanstd(mi_in_n_)/sqrt(length(mi_in_n_))
%
% nanmean(mi_out_n_)
% nanstd(mi_out_n_)/sqrt(length(mi_out_n_))
ind_neurons_resp=ind_it_neurons_ses(ix_h_both);
% save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\control_selection.mat','mi_in_n_','mi_out_n_','ind_neurons_resp')
% figure
% subplot(121)
% hold on
x_h_=[];x_h_=mi_out_n_;
y_h_=[];y_h_=mi_in_n_;
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
% plot([-1:.1:1],[-1:0.1:1])
p=signrank(x_h,y_h)
n=length(y_h);
aaaa=sum(x_h<y_h)/n
%cohen_d1=(nanmean(y_h)-nanmean(x_h))/sqrt((((length(y_h)-1)*var(y_h))+((length(x_h)-1)*var(x_h)))/(length(y_h)+length(x_h)-2))
cohen_d=((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h)
var_in_bts_test=squeeze(bootstrap_in(ix_h_both,:,:));
var_out_bts_test=squeeze(bootstrap_out(ix_h_both,:,:));

var_in_bts_test=var_in_bts_test(~isnan(x_h_)&~isnan(y_h_),:);
var_out_bts_test=var_out_bts_test(~isnan(x_h_)&~isnan(y_h_),:);

for bts=1:size(var_in_bts_test,1)
    
    try
        % pval_right(bts)=signrank(var_in_bts_test(bts,:),var_out_bts_test(bts,:));
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        pval_right(bts)=(nanmean(var_in_bts_test(bts,:))>prctile(var_out_bts_test(bts,:),95));
        pval_left(bts)=(nanmean(var_out_bts_test(bts,:))>prctile(var_in_bts_test(bts,:),95));
        
        % pval_left(bts)=signrank(var_out_bts_test(bts,:),var_in_bts_test(bts,:));
        
        
    catch
        pval_right(bts)=1;
        pval_left(bts)=1;
        
    end
end
% percent_of_increase=sum(pval_right<0.05)/length(pval_right)
% percent_of_decrease=sum(pval_left<0.05)/length(pval_right)
percent_of_increase=sum(pval_right==1)/length(pval_right)
percent_of_decrease=sum(pval_left==1)/length(pval_right)

eval(strcat('selectivity_mi_in_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

% scatter(x_h,y_h)
% scatter(x_h(m1_ind),y_h(m1_ind),'b*')
% scatter(x_h(m2_ind),y_h(m2_ind),'gs')
% export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5e',7)

% text( 0.2,0.3,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 0.2,0.2,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')
%
% set(gca,'fontsize',14,'fontweight','bold')
% xlabel('MI OUT')
% ylabel('MI IN')
%  axis([-0.4 0.6 -0.4 0.6])
% subplot(122)
% var_=[]; var_=-1*(y_h-x_h);
% % hist(var_,bin_num)
% bin_=0.05;
% nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
% [counts,centers]=hist(var_,nbin2);
% counts=counts/sum(counts);
% bar(centers,counts,1)
% line([0 0],ylim)
%
% xlim([-2 2.5])
%
% m=median(var_);
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m m],ylim,'linewidth',2,'color','r')
%
% m1=median(var_(m1_ind));
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m1 m1],ylim,'linewidth',2,'color','b')
%
% m2=median(var_(m2_ind));
% line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
% line([m2 m2],ylim,'linewidth',2,'color','g')
%
% if(p<0.05)
%     text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
%
%     text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
% end
% set(gca,'fontsize',14,'fontweight','bold')
% ylabel('% Pop.')
% xlim([-0.5 0.5])
%
% %%Fig behrad

%%Correlation acrsos all IT neurons IN

%  ix_h_all=[];ix_h_all=ind_it_neurons'&ix_h_tri_in&ix_h_tri_out&ind_pval;

% [m ~]=find(isnan(mi_high_out(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_low_out(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_high_in(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% [m ~]=find(isnan(mi_low_in(:,t_h>t1_delay&t_h<t2_delay)));m=unique(m);ix_h_all(m)=0;
% %

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    
    mi_all_n_=nanmean(mi_all_in(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_in(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_in(:,t_ind),2);
else
    mi_all_n_=mi_all_in_d;
    mi_low_n_=mi_low_in_d;
    mi_high_n_=mi_high_in_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
%  ix_h([9,178])=false;
%ix_h=ix_h_in;
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));
%
% figure('name','IN')
% subplot(121)
% title('IN')
% hold on
% scatter(x_h,y_h)
% scatter(x_h(m1_ind),y_h(m1_ind),'b*')
% scatter(x_h(m2_ind),y_h(m2_ind),'gs')
% export_file_to_excel([[ones(sum(m1_ind),1); 2*ones(sum(m2_ind),1)],x_h,y_h],'fig5f',7)

[r1, m ,b]=regression(x_h,y_h,'one');
[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

% plotregression(x_h,y_h)
%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coecloficients{2,1};
%%%lm.Coefficients{1,1};

xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

length(x_h)
% plot(xx,m*xx+b,'r')
% plot(xx,m_m1*xx+b_m1,'b')
% plot(xx,m_m2*xx+b_m2,'g')
% xlabel(' IT selectivity (a.u)')
% ylabel('Diff IT selectivity[High- Low](a.u)')
% ylim([-0.3 0.6])
%  xlim([0.1 1])
% ylim([-0.4 0.4])
% xlim([0.2 1])
% line([0.5 0.5],ylim)
% line(xlim,[0 0])
% set(gca,'fontsize',14,'fontweight','bold')
% text( 0.5,-0.05,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 0.5,-0.15,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')


%%% Correlation acrsos all It neurons OUt

mi_all_n_=[];mi_low_n_=[];mi_high_n_=[];
if ~interval
    mi_all_n_=nanmean(mi_all_out(:,t_ind),2);
    mi_low_n_=nanmean(mi_low_out(:,t_ind),2);
    mi_high_n_=nanmean(mi_high_out(:,t_ind),2);
else
    mi_all_n_=mi_all_out_d;
    mi_low_n_=mi_low_out_d;
    mi_high_n_=mi_high_out_d;
end
m1_ind=((ind_it_neurons_ses(ix_h_all))<52);
m2_ind=((ind_it_neurons_ses(ix_h_all))>51);
x_h_=[];x_h_=mi_all_n_(ix_h_all);
y_h_=[];y_h_=mi_high_n_(ix_h_all)-mi_low_n_(ix_h_all);

x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
m1_ind=m1_ind(~isnan(x_h_)&~isnan(y_h_));
m2_ind=m2_ind(~isnan(x_h_)&~isnan(y_h_));

[r, m ,b]=regression(x_h,y_h,'one');

[r_m1, m_m1 ,b_m1]=regression(x_h(m1_ind),y_h(m1_ind),'one');
[r_m2, m_m2 ,b_m2]=regression(x_h(m2_ind),y_h(m2_ind),'one');

%%%
y = y_h;
x = x_h;

tbl = table(x,y);
lm = fitlm(tbl,'linear');
% b=lm.Coefficients{1,1};
% m=lm.Coefficients{2,1};
xx=-0.3:0.1:1.5;
% [r,p]=corr(x_h,y_h,'type','kendal');
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,1001)
cohen_r2=(r^2)

length(x_h)


