function Fig1c_Firing_rate(resp)
%% In the name of Allah

%% Firing rate analysis
% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 1c, Fig S2a-d
% clear all
% close all
% clc

%% Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Data\gen_session_inf.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Data\gen_resp_pairs.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Data\gen_lfp_pairs.mat')
load('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\Color_codes')

win_=10;psth=@(x) SmoothN(1000*mean(x,1), win_);
ind_b=[201 550];
t_before_saccade=1600;t_after_saccade=100;

%% Setting

colors=[0 0 0;0.3 0.3 0.3;1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0; 0.5 0 0.5;0 0.7 0.7;1 1 0;0.7 0.7 0; 1 0 1;];
Color_monk=[0.9,0.9,0;0.5,0.5,0];

t_st_delay=1100;t_end_delay=1700;
N_sub=501;
N_bts=1001;
delay_si_selective_neurons=0;
si_selective=0;

%% Core calculator
cal_=1;
if cal_
    
    T_st=1;T_end=2100;
    ind_ses_fef=[];
    sui=1;
    for ss=1
                
        ss_p=ss;        
        n_su=size(resp(ss_p).FEF,2);
        for unini=1:n_su
            sig=[];  sig_=[];
            sig=resp(ss_p).FEF{unini}(:,T_st:T_end);
            N_trials=size(sig,1);
            for tt=1:size(sig,1)
                sig_(tt,:)=SmoothN(sig(tt,:)*1, 10);
            end
            
%             sig_m_=squeeze(mean(sig_(:,1:500),2));
%             sig_m=mean(sig_m_);sig_std=std(sig_m_);
%             sig_n=(sig_-sig_m)/sig_std;
%             
            for condi=1:2
                if condi<2; condi_h=[1,2,3];else condi_h=[1,2,3]+3; end
                
                ix_h_cr = []; ix_h_cr = logical(ismember(resp(ss_p).condition,condi_h))&(resp(ss_p).bhv==1)';
                ix_h_wr = []; ix_h_wr = logical(ismember(resp(ss_p).condition,condi_h))&(resp(ss_p).bhv==0)';
                
                N_trial_cr=sum(ix_h_cr);
                N_trial_wr=sum(ix_h_wr);
                
                var_cr_h = [];var_cr_h = sig_(ix_h_cr,:);
                var_wr_h = [];var_wr_h = sig_(ix_h_wr,:);
                
                
                mean_cr_h = mean(var_cr_h);                
                
                if N_trial_wr>1
                    mean_wr_h=mean(var_wr_h);
                end
                
                Val_cr_FEF_psth(sui,condi,:)=mean_cr_h;
                Val_wr_FEF_psth(sui,condi,:) = mean_wr_h;
                
                pval_da_fef(sui,condi)=signrank([mean(var_cr_h(:,ind_b(1):ind_b(2)),2);mean(var_wr_h(:,ind_b(1):ind_b(2)),2)],[mean(var_cr_h(:,t_st_delay:t_end_delay),2);mean(var_wr_h(:,t_st_delay:t_end_delay),2)],'tail','left');
                
            end
            
            var_h=[squeeze(Val_cr_FEF_psth(sui,:,:));squeeze(Val_wr_FEF_psth(sui,:,:))];
            var_b=mean2(var_h(:,ind_b(1):ind_b(2)));
            var_max=max(max(var_h));
            Val_cr_FEF_psth_norm(sui,:,:)= (squeeze(Val_cr_FEF_psth(sui,:,:))-var_b)/(var_max-var_b);
            Val_wr_FEF_psth_norm(sui,:,:)= (squeeze(Val_wr_FEF_psth(sui,:,:))-var_b)/(var_max-var_b);
            
            
            
            fr_in= mean2(sig_(logical(ismember(resp(ss_p).condition,[1 2 3])),t_st_delay:t_end_delay));
            fr_out= mean2(sig_(logical(ismember(resp(ss_p).condition,[4 5 6])),t_st_delay:t_end_delay));
            si=(fr_in-fr_out)/(fr_in+fr_out);
            for iteri=1:N_bts
                ind_h=randperm(N_trials,round(N_trials/2));
                fr_in_h= mean2(sig_(ind_h,t_st_delay:t_end_delay));
                ind_h=setdiff([1:N_trials],ind_h);
                fr_out_h= mean2(sig_(ind_h,t_st_delay:t_end_delay));
                si_h(iteri)=(fr_in_h-fr_out_h)/(fr_in_h+fr_out_h);
            end
            pval_dsi_fef(sui)=1-sum(si>si_h)/N_bts;
            
            sui=sui+1;
            ind_ses_fef=[ind_ses_fef ss];
            
        end
        
    end
    ind_ses_it=[];
    sui=1;
    for ss=1
              
        ss_p=ss;        
        n_su=size(resp(ss_p).IT,2);
        for unini=1:n_su
            sig=[];  sig_=[];
            sig=resp(ss_p).IT{unini}(:,T_st:T_end);
            N_trials=size(sig,1);
            for tt=1:size(sig,1)
                sig_(tt,:)=SmoothN(sig(tt,:)*1, 10);
            end
            
%             sig_m_=squeeze(mean(sig_(:,1:500),2));
%             sig_m=mean(sig_m_);sig_std=std(sig_m_);
%             sig_n=(sig_-sig_m)/sig_std;
%             
            obj_ix_h = [resp(ss_p).pref_obj_su(unini),resp(ss_p).npref_obj_su(unini)];
            obj_ix_h = [obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2)];
            obj_ix_h=[obj_ix_h];
            for condi=1:3
                condi_h=[obj_ix_h(condi) obj_ix_h(condi)+3];
                
                ix_h_cr = []; ix_h_cr = logical(ismember(resp(ss_p).condition,condi_h))&(resp(ss_p).bhv==1)';
                ix_h_wr = []; ix_h_wr = logical(ismember(resp(ss_p).condition,condi_h))&(resp(ss_p).bhv==0)';
                
                N_trial_cr=sum(ix_h_cr);
                N_trial_wr=sum(ix_h_wr);
                
                var_cr_h = [];var_cr_h = sig_(ix_h_cr,:);
                var_wr_h = [];var_wr_h = sig_(ix_h_wr,:);
                                                   
                mean_cr_h = mean(var_cr_h);
                               
                if N_trial_wr>1
                    mean_wr_h=mean(var_wr_h);
                end
                
                Val_cr_IT_psth(sui,condi,:)=mean_cr_h;
                Val_wr_IT_psth(sui,condi,:) = mean_wr_h;
                
                 pval_da_it(sui,condi)=signrank([mean(var_cr_h(:,ind_b(1):ind_b(2)),2);mean(var_wr_h(:,ind_b(1):ind_b(2)),2)],[mean(var_cr_h(:,t_st_delay:t_end_delay),2);mean(var_wr_h(:,t_st_delay:t_end_delay),2)],'tail','left');
                
%                  pval_da_it(sui,condi)=signrank([mean(var_cr_h(:,1:500),2);mean(var_wr_h(:,1:500),2)],[mean(var_cr_h(:,t1:t2),2);mean(var_wr_h(:,t1:t2),2)]);

                
            end
            
            var_h=[squeeze(Val_cr_IT_psth(sui,:,:));squeeze(Val_wr_IT_psth(sui,:,:));];
            var_b=mean2(var_h(:,ind_b(1):ind_b(2)));
            var_max=max(max(var_h));
            Val_cr_IT_psth_norm(sui,:,:)= (squeeze(Val_cr_IT_psth(sui,:,:))-var_b)/(var_max-var_b);
            Val_wr_IT_psth_norm(sui,:,:)= (squeeze(Val_wr_IT_psth(sui,:,:))-var_b)/(var_max-var_b);
            
            
            fr_in= mean2(sig_(logical(ismember(resp(ss_p).condition,[obj_ix_h(1) obj_ix_h(1)+3])),t_st_delay:t_end_delay));
            fr_out= mean2(sig_(logical(ismember(resp(ss_p).condition,[obj_ix_h(3) obj_ix_h(3)+3])),t_st_delay:t_end_delay));
            si=(fr_in-fr_out)/(fr_in+fr_out);
            for iteri=1:N_bts
                ind_h=randperm(N_trials,round(N_trials/3));
                fr_in_h= mean2(sig_(ind_h,t_st_delay:t_end_delay));
                
                ind_h_2=setdiff([1:N_trials],ind_h);
                ind_h_ind=randperm(length(ind_h_2),round(N_trials/3));
                ind_h=ind_h_2(ind_h_ind);
                fr_out_h= mean2(sig_(ind_h,t_st_delay:t_end_delay));
                si_h(iteri)=(fr_in_h-fr_out_h)/(fr_in_h+fr_out_h);
            end
            pval_dsi_it(sui)=1-sum(si>si_h)/N_bts;
            
            sui=sui+1;
            ind_ses_it=[ind_ses_it ss];
            
        end
        
    end
    
%         save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\Normalized_firing_rate_full.mat'...
%             ,'Val_cr_IT_psth_norm','Val_wr_IT_psth_norm'...
%             ,'Val_cr_FEF_psth_norm','Val_wr_FEF_psth_norm'...
%             ,'ind_ses_it','ind_ses_fef','pval_da_it','pval_da_fef','pval_dsi_it','pval_dsi_fef')
%     
%         save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_full.mat'...
%             ,'Val_cr_IT_psth','Val_wr_IT_psth'...
%             ,'Val_cr_FEF_psth','Val_wr_FEF_psth'...
%             ,'ind_ses_it','ind_ses_fef','pval_da_it','pval_da_fef','pval_dsi_it','pval_dsi_fef')
%     
%     
    
end
%% Load data

% load ('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\Firing Rate\Normalized_firing_rate.mat');
% load ('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\Firing Rate\firing_rate.mat');

psth_cr_it_n=Val_cr_IT_psth_norm;psth_wr_it_n=Val_wr_IT_psth_norm;
psth_cr_fef_n=Val_cr_FEF_psth_norm;psth_wr_fef_n=Val_wr_FEF_psth_norm;

psth_cr_it=Val_cr_IT_psth;psth_wr_it=Val_wr_IT_psth;
psth_cr_fef=Val_cr_FEF_psth;psth_wr_fef=Val_wr_FEF_psth;

t_h=[1:2100]-500;

%% Session selection for preview corect 
ind_ses=session_selection(1);
idx_it=find(ismember(ind_ses_it,ind_ses));
ind_ses_m1=setdiff(ind_ses,[52:92]);
ind_ses_m2=setdiff(ind_ses,[1:51;]);
it_m1=find(ismember(ind_ses_it,ind_ses_m1));
it_m2=find(ismember(ind_ses_it,ind_ses_m2));
ind_it_m1=(ind_ses_it<52);
ind_it_m2=(ind_ses_it>51);

ind_ses=session_selection(1);
idx_fef=find(ismember(ind_ses_fef,ind_ses));
ind_fef_m1=(ind_ses_fef<52);
ind_fef_m2=(ind_ses_fef>51);
ind_ses_m1=setdiff(ind_ses,[52:92;]);
ind_ses_m2=setdiff(ind_ses,[1:51;]);

fef_m1=find(ismember(ind_ses_fef,ind_ses_m1));
fef_m2=find(ismember(ind_ses_fef,ind_ses_m2));
%idx_fef=find(ismember(ind_ses_fef,ind_ses));

%
th_h=0.002;p_val_d=0.05;
ixx_it=idx_it;ixx_fef=idx_fef;

ih_s=ixx_it(find(abs(mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find(abs(mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ind_it_m1=ind_it_m1(ixx_it);
ind_it_m2=ind_it_m2(ixx_it);

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ind_fef_m1=ind_fef_m1(ixx_fef);
ind_fef_m2=ind_fef_m2(ixx_fef);

fef_m1=find(ismember(ind_ses_fef(ixx_fef),ind_ses_m1));
fef_m2=find(ismember(ind_ses_fef(ixx_fef),ind_ses_m2));

it_m1=find(ismember(ind_ses_it(ixx_it),ind_ses_m1));
it_m2=find(ismember(ind_ses_it(ixx_it),ind_ses_m2));

% 
% ixx_it_pa=ixx_it((pval_dsi_it(ixx_it)<p_val_d));
% ixx_fef_pa=ixx_fef((pval_dsi_fef(ixx_fef)<p_val_d));
%
if si_selective
ixx_it_pa=ixx_it((pval_dsi_it(ixx_it)<p_val_d));
ixx_fef_pa=ixx_fef((pval_dsi_fef(ixx_fef)<p_val_d));
else
ixx_it_pa=ixx_it((pval_da_it(ixx_it,1)<p_val_d));
ixx_fef_pa=ixx_fef((pval_da_fef(ixx_fef,1)<p_val_d));
% ixx_it_pa=setdiff(ixx_it_pa,[1:180]);
%      ixx_it_pa=setdiff(ixx_it_pa,[181:240]);
end


ixx_it_npa=setdiff(ixx_it,ixx_it_pa);
ixx_fef_npa=setdiff(ixx_fef,ixx_fef_pa);
%
if delay_si_selective_neurons
    ixx_it=ixx_it_pa;
    ixx_fef=ixx_fef_pa;
end

%% Fig. 1c Preview population activity
a=1;
figure('name','Fig. 1c, Rezayat.E, et al')
ax = subplot(2,4,1:4);    hold on;

for c=1:2;
    c1 = Color_loc(c,:);    
    niceplot(squeeze(psth_cr_fef_n(ixx_fef,c,:)),t_h/1000,win_,c1(1),c1(2),c1(3));
end
xlabel('Time from sample onset (sec.)');
ylabel('Normalized firing rate (a.u.)');
title(strcat('FEF'));
text(0.00,0.01,strcat(' N = ',num2str(length(ixx_fef))),'fontsize',14,'fontweight','bold');
line([0 0],ylim,'Color','k');line([0.300 0.300],ylim,'Color','k');line([1.300 1.300],ylim,'Color','k');
line([0.600 1.200],[-0.08 -0.08],'Color','k','Linewidth',4);
ax.XLim = [-0.300 1.600]; ax.YLim=[-0.08 0.5];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[0 0.5];
text(0.7,0.4,'In','color',Color_loc(1,:),'fontsize',20)
text(0.7,0.3,'Out','color',Color_loc(2,:),'fontsize',20)
 
ax = subplot(2,4,5:8);
    hold on;
for c=1:3;
    c1 = Color_obj(c,:);
   
    niceplot(squeeze(psth_cr_it_n(ixx_it,c,:)),t_h/1000,win_,c1(1),c1(2),c1(3));
    
    
end
xlabel('Time from sample onset (sec.)');
ylabel('Normalized firing rate (a.u)');
title(strcat('IT'));
text(0.00,0.01,strcat(' N = ',num2str(length(ixx_it))),'fontsize',14,'fontweight','bold');
line([0 0],ylim,'Color','k');line([0.300 0.300],ylim,'Color','k');line([1.300 1.300],ylim,'Color','k');
ax.XLim = [-0.300 1.600]; ax.YLim=[-0.08 0.5];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[0 0.5];
text(0.7,0.4,'Pref','color',Color_obj(1,:),'fontsize',20)
text(0.7,0.3,'Inter','color',Color_obj(2,:),'fontsize',20)
text(0.7,0.2,'NPref','color',Color_obj(3,:),'fontsize',20)

%% Stat for Selectivity
disp('%%% State for Selectivity %%%')
x_h=squeeze(mean(psth_cr_it_n(ixx_it,3,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT object selectivity, Pref vs. NPref, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF spatial selectivity, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

%% Stat compare to Baseline 
disp('%%% State for Delay activity  vs. Baseline %%%')
x_h=squeeze(mean(psth_cr_it_n(ixx_it,1,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Delay activity, Pref vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Delay activity, In vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_cr_it_n(ixx_it,3,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,3,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Delay activity, NPref vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Delay activity, out vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

%% Percent of neurons with delay activity
disp('%%% Percent of neurons with delay activity %%%')

disp (strcat('IT Delay activite neruons: ',' N',' =',num2str(length(ixx_it_pa)),' ; N all = ',...
    num2str(length(ixx_it)),'; Percent = ',num2str(length(ixx_it_pa)/length(ixx_it)*100),' % '));
disp (strcat('FEF Delay activite neruons: ',' N',' =',num2str(length(ixx_fef_pa)),' ; N all = ',...
    num2str(length(ixx_fef)),'; Percent = ',num2str(length(ixx_fef_pa)/length(ixx_fef)*100),' % '));

%% Session selection for comparing corect and Wrong
ind_ses=session_selection(6);
idx_it=find(ismember(ind_ses_it,ind_ses));
ind_ses_m1=setdiff(ind_ses,[52:92]);
ind_ses_m2=setdiff(ind_ses,[1:51]);
it_m1=find(ismember(ind_ses_it,ind_ses_m1));
it_m2=find(ismember(ind_ses_it,ind_ses_m2));

ind_ses=session_selection(5);
idx_fef=find(ismember(ind_ses_fef,ind_ses));
ind_ses_m1=setdiff(ind_ses,[52:92;]);
ind_ses_m2=setdiff(ind_ses,[1:51;]);
fef_m1=find(ismember(ind_ses_fef,ind_ses_m1));
fef_m2=find(ismember(ind_ses_fef,ind_ses_m2));
ind_fef_m1=(ind_ses_fef<52);
ind_fef_m2=(ind_ses_fef>51);

th_h=0.002;
ixx_it=idx_it;ixx_fef=idx_fef;
ih_s=ixx_it(find(abs(mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find(abs(mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ind_it_m1=(ind_ses_it(ixx_it)<52);
ind_it_m2=(ind_ses_it(ixx_it)>51);
it_m1=ixx_it(ind_it_m1);
it_m2=ixx_it(ind_it_m2);

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);

ind_fef_m1=(ind_ses_fef(ixx_fef)<52);
ind_fef_m2=(ind_ses_fef(ixx_fef)>51);
fef_m1=ixx_fef(ind_fef_m1);
fef_m2=ixx_fef(ind_fef_m2);

ixx_it_npa=setdiff(ixx_it,ixx_it_pa);
ixx_fef_npa=setdiff(ixx_fef,ixx_fef_pa);

%% Fig. 1c Preview population activity
figure('name','Fig. 1c, Rezayat.E, et al')
ax = subplot(2,4,1:4);    hold on;

c1 = Color_bhv(2,:);
niceplot(squeeze((psth_wr_fef_n(ixx_fef,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((psth_cr_fef_n(ixx_fef,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
xlabel('Time from sample onset (sec.)');
ylabel('Normalized firing rate (a.u.)');
title(strcat('FEF'));
text(0.600,0.2,strcat(' N = ',num2str(length(ixx_fef))),'fontsize',14,'fontweight','bold');
line([0 0],ylim,'Color','k');line([0.300 0.300],ylim,'Color','k');line([1.300 1.300],ylim,'Color','k');
ax.XLim = [-0.300 1.600]; ax.YLim=[-0.08 0.5];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[0 0.5];
text(0.7,0.4,'Cr','color',Color_bhv(1,:),'fontsize',20)
text(0.7,0.3,'Wr','color',Color_bhv(2,:),'fontsize',20)
ax=subplot(2,4,5:8);
hold on
c1 = Color_bhv(2,:);
niceplot(squeeze((psth_wr_it_n(ixx_it,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((psth_cr_it_n(ixx_it,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
xlabel('Time from sample onset (sec.)');
ylabel('Normalized firing rate (a.u.)');
title(strcat('IT'));
text(0.600,0.2,strcat(' N = ',num2str(length(ixx_it))),'fontsize',14,'fontweight','bold');
line([0 0],ylim,'Color','k');line([0.300 0.300],ylim,'Color','k');line([1.300 1.300],ylim,'Color','k');
ax.XLim = [-0.300 1.600]; ax.YLim=[-0.08 0.5];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[0 0.5];

%% Fig. S2a Scatter Firing rate Cr and Wr FEF
y_h= squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
x_h= squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));
ind_m1_=ind_fef_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_fef_m2((~isnan(x_h_)&~isnan(y_h_)));

figure('name','Fig. S2b, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
scatter(x_h(ind_m1_),y_h(ind_m1_),'b*')
scatter(x_h(ind_m2_),y_h(ind_m2_),'gs')
p=signrank(x_h,y_h);
n=length(x_h);
if p>0.00099
    text( 0.20,0.30,strcat('p = ',num2str(p,1)),'fontsize',14,'fontweight','bold')
    
else
    text( 0.20,0.30,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
end
text( 0.20,0.20,strcat('n = ',num2str(n,3)),'fontsize',14,'fontweight','bold')

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Norm. Firing rate (a.u.) [Wr]');
ylabel('Norm. Firing rate (a.u.) [Cr]');
axis([-0.1 0.6 -0.1 0.6])
ax.XTick=[-0.1 0.6];
ax.YTick=[-0.1 0.6];

%% State FEF Firing rate  Cr vs. Wr 
try
disp('%%% State FEF Firing rate  Cr vs. Wr  %%%')

x_h=squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_fef_n(ixx_fef(ind_fef_m1),1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef(ind_fef_m1),1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_fef_n(ixx_fef(ind_fef_m2),1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef(ind_fef_m2),1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end

%% Fig. S2b Scatter Cr and Wr IT
y_h= squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
x_h= squeeze(mean(psth_wr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));

ind_m1_=ind_it_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_it_m2((~isnan(x_h_)&~isnan(y_h_)));

p=signrank(x_h,y_h);
n=length(x_h);
 m_it_cr_wr_fr =  [nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h) ,n, p];

figure('name','Fig. S2b, Rezayat.E, et al')
hold on
scatter(x_h(ind_m1_),y_h(ind_m1_),'b*')
scatter(x_h(ind_m2_),y_h(ind_m2_),'gs')
if p>0.00099
    text( 0.20,0.30,strcat('p = ',num2str(p,3)),'fontsize',14,'fontweight','bold')
    
else
    text( 0.20,0.30,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
end
text( 0.20,0.20,strcat('n = ',num2str(n,3)),'fontsize',14,'fontweight','bold')

plot([-10:0.1:150],[-10:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')

axis([-0.2 0.6 -0.2 0.6]) 

xlabel('Norm. Firing rate (a.u.) [Wr]');
ylabel('Norm. Firing rate (a.u.) [Cr]');
ax.XTick=[-0.2 0.6];
ax.YTick=[-0.2 0.6];

%% State IT Firing rate  Cr vs. Wr 
try
disp('%%% State IT Firing rate  Cr vs. Wr  %%%')
x_h=squeeze(mean(psth_wr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_it_n(ixx_it(ind_it_m1),1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it(ind_it_m1),1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_it_n(ixx_it(ind_it_m2),1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it(ind_it_m2),1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end

%% Fig. S2c Scatter FEF Spatial Selectivity Index Cr and Wr 
x_in=squeeze(mean(psth_cr_fef(ixx_fef,1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_cr_fef(ixx_fef,2,t_st_delay:t_end_delay),3));
y_h = [];y_h=(x_in-x_out)./(x_in+x_out);

x_in=squeeze(mean(psth_wr_fef(ixx_fef,1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_wr_fef(ixx_fef,2,t_st_delay:t_end_delay),3));
x_h = [];x_h=(x_in-x_out)./(x_in+x_out);

x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));
ind_m1_=ind_fef_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_fef_m2((~isnan(x_h_)&~isnan(y_h_)));

figure('name','Fig. S2c, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
scatter(x_h(ind_m1_),y_h(ind_m1_),'b*')
scatter(x_h(ind_m2_),y_h(ind_m2_),'gs')
p=signrank(x_h,y_h);
n=length(x_h);
if p>0.00099
    text( 0.20,0.30,strcat('p = ',num2str(p,1)),'fontsize',14,'fontweight','bold')
    
else
    text( 0.20,0.30,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
end
text( 0.20,0.20,strcat('n = ',num2str(n,3)),'fontsize',14,'fontweight','bold')

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Selectivity index (a.u.) [Wr]');
ylabel('Selectivity index (a.u.) [Cr]');
axis([-0.5 1 -0.5 1])
ax.XTick=[-0.5 1];
ax.YTick=[-0.5 1];

%% State for FEF Saptial Selectivity Index  Cr vs. Wr
try
disp('%%% State for FEF Saptial Selectivity Index  Cr vs. Wr  %%%')

x_in=squeeze(mean(psth_cr_fef(ixx_fef,1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_cr_fef(ixx_fef,2,t_st_delay:t_end_delay),3));
y_h = [];y_h=(x_in-x_out)./(x_in+x_out);

x_in=squeeze(mean(psth_wr_fef(ixx_fef,1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_wr_fef(ixx_fef,2,t_st_delay:t_end_delay),3));
x_h = [];x_h=(x_in-x_out)./(x_in+x_out);

p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Saptial Selectivity Index  Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_in=squeeze(mean(psth_cr_fef(ixx_fef(ind_fef_m1),1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_cr_fef(ixx_fef(ind_fef_m1),2,t_st_delay:t_end_delay),3));
y_h = [];y_h=(x_in-x_out)./(x_in+x_out);

x_in=squeeze(mean(psth_wr_fef(ixx_fef(ind_fef_m1),1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_wr_fef(ixx_fef(ind_fef_m1),2,t_st_delay:t_end_delay),3));
x_h = [];x_h=(x_in-x_out)./(x_in+x_out);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Saptial Selectivity Index  Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_in=squeeze(mean(psth_cr_fef(ixx_fef(ind_fef_m2),1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_cr_fef(ixx_fef(ind_fef_m2),2,t_st_delay:t_end_delay),3));
y_h = [];y_h=(x_in-x_out)./(x_in+x_out);

x_in=squeeze(mean(psth_wr_fef(ixx_fef(ind_fef_m2),1,t_st_delay:t_end_delay),3));
x_out=squeeze(mean(psth_wr_fef(ixx_fef(ind_fef_m2),2,t_st_delay:t_end_delay),3));
x_h = [];x_h=(x_in-x_out)./(x_in+x_out);

p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Saptial Selectivity Index  Cr vs. Wrr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Fig. S2d Scatter IT Object Selectivity Index Cr and Wr 
x_pref=squeeze(mean(psth_cr_it(ixx_it,1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_cr_it(ixx_it,3,t_st_delay:t_end_delay),3));
y_h=(x_pref-x_npref)./(x_pref+x_npref);

x_pref=squeeze(mean(psth_wr_it(ixx_it,1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_wr_it(ixx_it,3,t_st_delay:t_end_delay),3));
x_h=(x_pref-x_npref)./(x_pref+x_npref);

x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));
ind_m1_=ind_it_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_it_m2((~isnan(x_h_)&~isnan(y_h_)));

figure('name','Fig. S2c, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
scatter(x_h(ind_m1_),y_h(ind_m1_),'b*')
scatter(x_h(ind_m2_),y_h(ind_m2_),'gs')
p=signrank(x_h,y_h);
n=length(x_h);
if p>0.00099
    text( 0.20,0.30,strcat('p = ',num2str(p,1)),'fontsize',14,'fontweight','bold')
    
else
    text( 0.20,0.30,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
end
text( 0.20,0.20,strcat('n = ',num2str(n,3)),'fontsize',14,'fontweight','bold')

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Selectivity index (a.u.) [Wr]');
ylabel('Selectivity index (a.u.) [Cr]');
axis([-0.5 1 -0.5 1])
ax.XTick=[-0.5 1];
ax.YTick=[-0.5 1];

%% State for IT  Object  Selectivity Index  Cr vs. Wr
try
disp('%%% State for IT Object Selectivity Index  Cr vs. Wr  %%%')

x_pref=squeeze(mean(psth_cr_it(ixx_it,1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_cr_it(ixx_it,3,t_st_delay:t_end_delay),3));
y_h=(x_pref-x_npref)./(x_pref+x_npref);

x_pref=squeeze(mean(psth_wr_it(ixx_it,1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_wr_it(ixx_it,3,t_st_delay:t_end_delay),3));
x_h=(x_pref-x_npref)./(x_pref+x_npref);

p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Object Selectivity Index  Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_pref=squeeze(mean(psth_cr_it(ixx_it(ind_it_m1),1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_cr_it(ixx_it(ind_it_m1),3,t_st_delay:t_end_delay),3));
y_h=(x_pref-x_npref)./(x_pref+x_npref);

x_pref=squeeze(mean(psth_wr_it(ixx_it(ind_it_m1),1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_wr_it(ixx_it(ind_it_m1),3,t_st_delay:t_end_delay),3));
x_h=(x_pref-x_npref)./(x_pref+x_npref);

p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Object Selectivity Index  Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_pref=squeeze(mean(psth_cr_it(ixx_it(ind_it_m2),1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_cr_it(ixx_it(ind_it_m2),3,t_st_delay:t_end_delay),3));
y_h=(x_pref-x_npref)./(x_pref+x_npref);

x_pref=squeeze(mean(psth_wr_it(ixx_it(ind_it_m2),1,t_st_delay:t_end_delay),3));
x_npref=squeeze(mean(psth_wr_it(ixx_it(ind_it_m2),3,t_st_delay:t_end_delay),3));
x_h=(x_pref-x_npref)./(x_pref+x_npref);

p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Saptial Selectivity Index  Cr vs. Wrr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Stat for Target Period Fring rate Cr vs. Wr  
try
t_st_sac=1450;t_end_sac=2100;
disp('%%% Stat for Target Period Fring rate Cr vs. Wr   %%%')
x_h=squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_sac:t_end_sac),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_sac:t_end_sac),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, Target period; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_it_n(ixx_it,1,t_st_sac:t_end_sac),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_sac:t_end_sac),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, Target period; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Stat for visual Period Fring rate Cr vs. Wr  
try
t_st_sac=550;t_end_sac=850;
disp('%%% Stat for visual Period Fring rate Cr vs. Wr   %%%')
x_h=squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_sac:t_end_sac),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_sac:t_end_sac),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, visual period; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(mean(psth_wr_it_n(ixx_it,1,t_st_sac:t_end_sac),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_sac:t_end_sac),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, visual period; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

catch
end