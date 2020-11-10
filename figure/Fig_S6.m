%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. S62,  Control for Location coding in FEF
clear all
close all
clc

%% Settings
ind_b=[1 500];
Ns=92;
win_=300;step=10;
Shuffle=0;
N_shuff=30;
inf_cal=@(x,y) ROC_(x,y);

%% Rate categorized calculate

load([cd,'\Mat\PPL\ROC_FEF_interval.mat'])

%% Load PPL 

% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\PPL_full_shuff_1001-shuff.mat'));
load(strcat(cd,'\Mat\PPL\PPL_full_shuff_1001.mat'));

freqs2use_plv = freqs2use;
NS = size(plv_cr,1);
shuffle_corrected = 1;
if shuffle_corrected
    plv_cr = plv_cr_n;plv_wr = plv_wr_n;
end

%%Baseline Correction
baseline_corrected = 1;
ind_b = [41 110];
if baseline_corrected
    plv_cr_ = plv_cr;plv_wr_ = plv_wr;
    val_b = nanmean(plv_cr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    plv_cr = (plv_cr_-val_b)./val_sd;
    
    val_b = nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    plv_wr = (plv_wr_-val_b)./val_sd;
end


t = 1:2100;
t = decimate(t,5,'fir');
t_h = t-500;

%%
t1=600;t2=1200;
freq_bands = [3.9,7.9,21.9,35.9,50,70;8.1,15.1,28.1,50.1,130.1,130.1];
 
f=freqs2use_plv;
 f1=freq_bands(1,3);  f2=freq_bands(2,3);
 ind_ses=session_selection(2);
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                        % Sessions from 52-86 is for Monkey 2
  ind_ses=ind_ses_m1;

    var_cr1=squeeze(plv_cr(ind_ses,1,:,:))-squeeze(plv_wr(ind_ses,1,:,:));   
    y_0=(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,ind_b(1):ind_b(2)),3),2)));
    y_h=(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    

x_h=mi_all_mu((ind_ses),1);
ihs=(x_h>=0.5);
x_h=x_h(ihs);
 y_h=y_h(ihs);
ind_h_m2=ismember(ind_ses(ihs),[52:86]); 
ind_h_m1=ismember(ind_ses(ihs),[1:51]); 
p=signrank(x_h,y_h);
n=length(x_h);

%% control for location selectivity and PPL
ind_ses=session_selection(4);

ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                        % Sessions from 52-86 is for Monkey 2
%   ind_ses=ind_ses_m2;
y_h=[];
var_cr1=squeeze(plv_cr(:,1,:,:))-squeeze(plv_wr(:,1,:,:));
y_h=(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));

x_h=mi_all_both(ismember(ind_ses_fef,ind_ses),1);
ihs=(x_h>=-50.5);
ind_h=ind_ses_fef(ismember(ind_ses_fef,ind_ses));
x_h=x_h(ihs);
y_h_=y_h(ind_h(ihs));y_h=y_h_;

ind_h_m2=ismember(ind_h(ihs),[52:92]); 
ind_h_m1=ismember(ind_h(ihs),[1:51]); 
n=length(x_h);

figure('name','Fig S6 Rezayat.E, et al')

hold on
scatter(x_h(ind_h_m1),y_h(ind_h_m1),'b*')
scatter(x_h(ind_h_m2),y_h(ind_h_m2),'gs')
export_file_to_excel([[ones(sum(ind_h_m1),1); 2*ones(sum(ind_h_m2),1)],x_h,y_h],'figS6',11)

% plot([-100:250],[-100:250],'k')
ylabel('Diff PPL [Cr-Wr] (a.u.)')
xlabel('AUC FEF location sel. (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
% axis([0.1 0.25 0.1 0.25])
[r p]= correlation(x_h,y_h,1000);
[r, m ,b]=regression(x_h,y_h,'one');
% line([0.5 0.5],ylim)
n=length(x_h);
xx=0:0.1:1;
plot(xx,m*xx+b,'k')
disp (strcat(' All; ',' r=',num2str(r),' =',...
    '; n = ',num2str(n),' p = ',num2str(p)));
text( 0.2,2.5,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.2,2.0,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.2,1.5,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')

[r p]= correlation(x_h(ind_h_m1),y_h(ind_h_m1),1000);
[r, m ,b]=regression(x_h(ind_h_m1),y_h(ind_h_m1),'one');
line([0.5 0.5],ylim)
n=length(x_h(ind_h_m1));
xx=0:0.1:1;
plot(xx,m*xx+b,'b')


disp (strcat(' M1; ',' r=',num2str(r),' =',...
    '; n = ',num2str(n),' p = ',num2str(p)));


[r p]= correlation(x_h(ind_h_m2),y_h(ind_h_m2),1000);
[r, m ,b]=regression(x_h(ind_h_m2),y_h(ind_h_m2),'one');
line([0.5 0.5],ylim)
n=length(x_h(ind_h_m2));
xx=0:0.1:1;
plot(xx,m*xx+b,'g')


disp (strcat(' M2; ',' r=',num2str(r),' =',...
    '; n = ',num2str(n),' p = ',num2str(p)));


