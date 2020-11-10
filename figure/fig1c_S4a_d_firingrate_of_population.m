%% In the name of Allah

%% Firing rate analysis
% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 1c, Fig S2a-d
clear all
close all
clc

%% Load data
load([ cd ,'\Data\gen_session_inf.mat'])
load([ cd ,'\Data\gen_resp_pairs.mat'])
load([ cd ,'\Data\gen_lfp_pairs.mat'])
load([ cd ,'\Data\Color_codes'])

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

%% Load data

load ([ cd ,'\Mat\Firing Rate\Normalized_firing_rate.mat']);
load ([ cd ,'\Mat\Firing Rate\firing_rate.mat']);

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
 x_export=t_h';

for c=1:2;
    c1 = Color_loc(c,:);    
    niceplot(squeeze(psth_cr_fef_n(ixx_fef,c,:)),t_h/1000,win_,c1(1),c1(2),c1(3));
%         export_file_to_excel([t_h/1000;squeeze(psth_cr_fef_n(ixx_fef,c,:)) ],['fig1cfef' num2str(c)] )
x_export=[x_export, nanmean(squeeze(psth_cr_fef_n(ixx_fef,c,:)))'];
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
    
%     export_file_to_excel([t_h/1000;squeeze(psth_cr_it_n(ixx_it,c,:)) ],['fig1cit' num2str(c)] )
x_export=[x_export, nanmean(squeeze(psth_cr_it_n(ixx_it,c,:)))'];

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
 disp (strcat('IT object selectivity, Pref vs. NPref,All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF spatial selectivity, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
 disp (strcat('FEF spatial selectivity, In vs. Out, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

%% Stat compare to Baseline 
disp('%%% State for Delay activity  vs. Baseline %%%')
x_h=squeeze(mean(psth_cr_it_n(ixx_it,1,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Delay activity, Pref vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
 disp (strcat('IT Delay activity, Pref vs. baseline, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Delay activity, In vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
disp (strcat('FEF Delay activity, In vs. baseline, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

x_h=squeeze(mean(psth_cr_it_n(ixx_it,3,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,3,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Delay activity, NPref vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
disp (strcat('IT Delay activity, NPref vs. baseline, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

x_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,ind_b(1):ind_b(2)),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,2,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF Delay activity, out vs. baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
disp (strcat('FEF Delay activity, out vs. baseline, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

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
%     export_file_to_excel([t_h/1000;squeeze(psth_cr_fef_n(ixx_fef,1,:)) ],['fig1cfefcr' ] )
x_export=[x_export, nanmean(squeeze(psth_wr_fef_n(ixx_fef,1,:)))'];

c1 = Color_bhv(1,:);
niceplot(squeeze((psth_cr_fef_n(ixx_fef,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
%     export_file_to_excel([t_h/1000;squeeze(psth_wr_fef_n(ixx_fef,1,:)) ],['fig1cfefwr' ] )
x_export=[x_export, nanmean(squeeze(psth_wr_it_n(ixx_it,1,:)))'];

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
%     export_file_to_excel([t_h/1000;squeeze(psth_cr_it_n(ixx_it,1,:)) ],['fig1citcr' ] )

c1 = Color_bhv(1,:);
niceplot(squeeze((psth_cr_it_n(ixx_it,1,:))),t_h/1000,win_,c1(1),c1(2),c1(3))
%     export_file_to_excel([t_h/1000;squeeze(psth_wr_it_n(ixx_it,1,:)) ],['fig1citwr' ] )

xlabel('Time from sample onset (sec.)');
ylabel('Normalized firing rate (a.u.)');
title(strcat('IT'));
text(0.600,0.2,strcat(' N = ',num2str(length(ixx_it))),'fontsize',14,'fontweight','bold');
line([0 0],ylim,'Color','k');line([0.300 0.300],ylim,'Color','k');line([1.300 1.300],ylim,'Color','k');
ax.XLim = [-0.300 1.600]; ax.YLim=[-0.08 0.5];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[0 0.5];
% 
export_file_to_excel(x_export,'fig1c',5)

%% Fig. S4a Scatter Firing rate Cr and Wr FEF
y_h= squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
x_h= squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));
ind_m1_=ind_fef_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_fef_m2((~isnan(x_h_)&~isnan(y_h_)));

figure('name','Fig. S4a, Rezayat.E, et al')
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
%   export_file_to_excel([x_h y_h],['figS4a' ] )
export_file_to_excel([[ones(sum(ind_m1_),1); 2*ones(sum(ind_m2_),1)],x_h,y_h],'figS4a',1)

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Norm. Firing rate (a.u.) [Wr]');
ylabel('Norm. Firing rate (a.u.) [Cr]');
axis([-0.1 0.6 -0.1 0.6])
ax.XTick=[-0.1 0.6];
ax.YTick=[-0.1 0.6];

%% State FEF Firing rate  Cr vs. Wr 
disp('%%% State FEF Firing rate  Cr vs. Wr  %%%')

x_h=squeeze(mean(psth_wr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_fef_n(ixx_fef,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
disp (strcat('FEF In, Cr vs. Wr, All, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

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

%% Fig. S4b Scatter Cr and Wr IT
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

figure('name','Fig. S4b, Rezayat.E, et al')
hold on
scatter(x_h(ind_m1_),y_h(ind_m1_),'b*')
scatter(x_h(ind_m2_),y_h(ind_m2_),'gs')
if p>0.00099
    text( 0.20,0.30,strcat('p = ',num2str(p,3)),'fontsize',14,'fontweight','bold')
    
else
    text( 0.20,0.30,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
end
text( 0.20,0.20,strcat('n = ',num2str(n,3)),'fontsize',14,'fontweight','bold')
%   export_file_to_excel([x_h y_h],['figS4b' ] )
export_file_to_excel([[ones(sum(ind_m1_),1); 2*ones(sum(ind_m2_),1)],x_h,y_h],'figS4b',1)

plot([-10:0.1:150],[-10:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')

axis([-0.2 0.6 -0.2 0.6]) 

xlabel('Norm. Firing rate (a.u.) [Wr]');
ylabel('Norm. Firing rate (a.u.) [Cr]');
ax.XTick=[-0.2 0.6];
ax.YTick=[-0.2 0.6];

%% State IT Firing rate  Cr vs. Wr 
disp('%%% State IT Firing rate  Cr vs. Wr  %%%')
x_h=squeeze(mean(psth_wr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
y_h=squeeze(mean(psth_cr_it_n(ixx_it,1,t_st_delay:t_end_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
disp (strcat('IT Pref, Cr vs. Wr, All;, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

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

%% Fig. S4c Scatter FEF Spatial Selectivity Index Cr and Wr 
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

figure('name','Fig. S4c, Rezayat.E, et al')
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
%   export_file_to_excel([x_h y_h],['figS4c' ] )
export_file_to_excel([[ones(sum(ind_m1_),1); 2*ones(sum(ind_m2_),1)],x_h,y_h],'figS4c',1)

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Selectivity index (a.u.) [Wr]');
ylabel('Selectivity index (a.u.) [Cr]');
axis([-0.5 1 -0.5 1])
ax.XTick=[-0.5 1];
ax.YTick=[-0.5 1];

%% State for FEF Saptial Selectivity Index  Cr vs. Wr
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
disp (strcat('FEF Saptial Selectivity Index  Cr vs. Wr All;, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

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

%% Fig. S4d Scatter IT Object Selectivity Index Cr and Wr 
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

figure('name','Fig. S4d, Rezayat.E, et al')
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
%   export_file_to_excel([x_h y_h],['figS4d' ] )
export_file_to_excel([[ones(sum(ind_m1_),1); 2*ones(sum(ind_m2_),1)],x_h,y_h],'figS4d',1)

plot([-1:0.1:150],[-1:0.1:150],'k:')
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Selectivity index (a.u.) [Wr]');
ylabel('Selectivity index (a.u.) [Cr]');
axis([-0.5 1 -0.5 1])
ax.XTick=[-0.5 1];
ax.YTick=[-0.5 1];

%% State for IT  Object  Selectivity Index  Cr vs. Wr
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
disp (strcat('IT Object Selectivity Index  Cr vs. Wr, All;, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

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

%% Stat for Target Period Fring rate Cr vs. Wr  
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

%% Stat for visual Period Fring rate Cr vs. Wr  
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

