%% Set Paths
clear all
clc
close all

%% Settings
M1=0;
M2=0;

normalize_= @(x) x/nansum(x);
plv=@(x) abs(nanmean(exp(1j*x)));
control=0;
Color_loc=[1,0,0;0,0,1];
Color_obj=[1,0.5,0;0.5,0.5,0.5;0,1,0.5];
Color_bhv=[1,1,0;0,1,1];
win_=1;
t1=1100;t2=1700;N_bts=1001;

freq_bands=[3.9,7.9,21.9,34.9,49.9,49.9;8.1,15.1,28.1,50.1,130.1,150.1];
freqs2use=4:60;

%% Load
load([cd,'\Mat\SPL\spl_full_withinarea_pairs.mat']);

load([cd,'\Mat\SPL\spl_full_intra_area_pairs.mat']);

%% Session selection for comparing corect and Wrong
load ([cd '\Mat\Firing Rate\firing_rate_full.mat']);

% load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_500msbaseline_alldata.mat');

psth_cr_it=Val_cr_IT_psth;psth_wr_it=Val_wr_IT_psth;
psth_cr_fef=Val_cr_FEF_psth;psth_wr_fef=Val_wr_FEF_psth;
%%
%status = 1 All , status = 2 SFN; status = 3 performance
ind_ses_itphase=session_select_spl_after_review(6);
%  ind_ses_itphase=setdiff(ind_ses_itphase,[74:92]);
if M1
    ind_ses_itphase=setdiff(ind_ses_itphase,[52:92]);
elseif M2
    ind_ses_itphase=setdiff(ind_ses_itphase,[1:51]);
    
end
load multiunit_preferance_pval.mat

p=0.05;
ind_nonslective=find((p_rate>p)&(p_baseline>p));

idx_it_spl=find(ismember(ind_ses_it2,ind_ses_itphase));

m1_it_ind=((ind_ses_it2)<52);
m2_it_ind=((ind_ses_it2)>51);
idx_fef_spl=find(ismember(ind_ses_fef2,ind_ses_itphase));
m1_fef_ind=((ind_ses_fef2)<52);
m2_fef_ind=((ind_ses_fef2)>51);

idx_fef=idx_fef_spl;
idx_it=idx_it_spl;

%
th_h=0.002;p_val_d=0.05;
ixx_it=idx_it;ixx_fef=idx_fef;


ixx_fef=idx_fef;

ih_s=ixx_it(find((mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find((mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it_itphase=setdiff(ixx_it,ih_s);


%%%% Control for selecting one unit from each session
if control
    idx_control=zeros(length(ind_ses_it2),1);
    ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_itphase)=1;
    
    for ss=1:92
        ind_h= find(ismember(ind_ses_it2,ss)'&ind_resp);
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    ixx_it_itphase=find(ind_resp&idx_control);
end
it_itphase_ind_m1=(ind_ses_it2(ixx_it_itphase)<52);
it_itphase_ind_m2=(ind_ses_it2(ixx_it_itphase)>51);


ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef_itphase=setdiff(ixx_fef,ih_s);
fef_itphase_ind_m1=(ind_ses_fef2(ixx_fef_itphase)<52);
fef_itphase_ind_m2=(ind_ses_fef2(ixx_fef_itphase)>51);
%%
ind_ses_fefphase=session_select_spl_after_review(8);   %% FOr preferance put it 7

if M1
    ind_ses_fefphase=setdiff(ind_ses_fefphase,[52:92]);
elseif M2
    ind_ses_fefphase=setdiff(ind_ses_fefphase,[1:51]);
    
end

idx_it_spl=find(ismember(ind_ses_it2,ind_ses_fefphase));
idx_fef_spl=find(ismember(ind_ses_fef2,ind_ses_fefphase));

idx_fef=idx_fef_spl;
idx_it=idx_it_spl;

%
th_h=0.002;p_val_d=0.05;
ixx_it=idx_it;ixx_fef=idx_fef;


ixx_fef=idx_fef;

ih_s=ixx_it(find((mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find((mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it_fefphase=setdiff(ixx_it,ih_s);

%%%% Control for selecting one unit from each session
if control
    idx_control=zeros(length(ind_ses_it2),1);
    ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_fefphase)=1;
    
    for ss=1:92
        ind_h= find(ismember(ind_ses_it2,ss)'&ind_resp);
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    ixx_it_fefphase=find(ind_resp&idx_control);
end
it_fefphase_ind_m1=(ind_ses_it2(ixx_it_fefphase)<52);
it_fefphase_ind_m2=(ind_ses_it2(ixx_it_fefphase)>51);

%

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef_fefphase=setdiff(ixx_fef,ih_s);
fef_fefphase_ind_m1=(ind_ses_fef2(ixx_fef_fefphase)<52);
fef_fefphase_ind_m2=(ind_ses_fef2(ixx_fef_fefphase)>51);

%
noramlized=0;
Color_loc=[1,0,0;0,0,1];
Color_obj=[0,0,0;0.5,0,0.5;0,1,0];
Color_bhv=[0.7,0.3,0;0,0.7,0.7];
%%
%status = 1 All , status = 2 SFN; status = 3 performance
ind_ses_itphase=session_select_spl_after_review(2);
%  ind_ses_itphase=setdiff(ind_ses_itphase,[74:92]);
if M1
    ind_ses_itphase=setdiff(ind_ses_itphase,[52:92]);
elseif M2
    ind_ses_itphase=setdiff(ind_ses_itphase,[1:51]);
    
end
load multiunit_preferance_pval.mat

p=0.05;
ind_nonslective=find((p_rate>p)&(p_baseline>p));

%   ind_ses_itphase=setdiff(ind_ses_itphase,ind_nonslective);

% ind_ses_it_rate=ind_ses_it(ismember(ind_ses_it,ind_ses_it2));
% ind_ses_fef_rate=ind_ses_fef(ismember(ind_ses_fef,ind_ses_fef2));
idx_it_spl=find(ismember(ind_ses_it2,ind_ses_itphase));

m1_it_ind=((ind_ses_it2)<52);
m2_it_ind=((ind_ses_it2)>51);
idx_fef_spl=find(ismember(ind_ses_fef2,ind_ses_itphase));
m1_fef_ind=((ind_ses_fef2)<52);
m2_fef_ind=((ind_ses_fef2)>51);
% idx_it_rate=find(ismember(ind_ses_it_rate,ind_ses_itphase));
% idx_fef_rate=find(ismember(ind_ses_fef_rate,ind_ses_itphase));
% idx_fef=find(ismember(idx_fef_spl,idx_fef_rate));
% idx_it=find(ismember(idx_it_spl,idx_it_rate));
%idx_fef=find(ismember(ind_ses_fef,ind_ses));
idx_fef=idx_fef_spl;
idx_it=idx_it_spl;

%
th_h=0.002;p_val_d=0.05;
ixx_it=idx_it;ixx_fef=idx_fef;


% ixx_it_=find(p_val_it(:,1)<0.05);
%  ixx_it=find(pval_dsi_it(idx_it_rate)<0.05);

% th_h=0.002;p_val_d=0.05;
%   ixx_it=idx_it;
ixx_fef=idx_fef;

ih_s=ixx_it(find((mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find((mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it_itphase=setdiff(ixx_it,ih_s);


%%%% Control for selecting one unit from each session 
if control
idx_control=zeros(length(ind_ses_it2),1);
ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_itphase)=1;

for ss=1:92
  ind_h= find(ismember(ind_ses_it2,ss)'&ind_resp);
  if ~isempty(ind_h)
  idx_control(ind_h(randi(length(ind_h),1)))=1;
  end
end
ixx_it_itphase=find(ind_resp&idx_control);
end
it_itphase_ind_m1=(ind_ses_it2(ixx_it_itphase)<52);
it_itphase_ind_m2=(ind_ses_it2(ixx_it_itphase)>51);

% 
ind_ses_fefphase=session_select_spl_after_review(4);   %% FOr preferance put it 4
%  ind_ses_fefphase=setdiff(ind_ses_fefphase,[74:92]);
%  ind_ses=setdiff(ind_ses,[52:73]);
% %     ind_ses_fefphase=setdiff(ind_ses_fefphase,[1:51]);
%    ind_ses_fefphase=setdiff(ind_ses_fefphase,[52:92]);
if M1
    ind_ses_fefphase=setdiff(ind_ses_fefphase,[52:92]);
elseif M2
    ind_ses_fefphase=setdiff(ind_ses_fefphase,[1:51]);
    
end

% ind_ses_it_rate=ind_ses_it(ismember(ind_ses_it,ind_ses_it2));
% ind_ses_fef_rate=ind_ses_fef(ismember(ind_ses_fef,ind_ses_fef2));
idx_fef_spl=find(ismember(ind_ses_fef2,ind_ses_fefphase));

idx_fef=idx_fef_spl;
%
th_h=0.002;p_val_d=0.05;
ixx_fef=idx_fef;

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef_fefphase=setdiff(ixx_fef,ih_s);
fef_fefphase_ind_m1=(ind_ses_fef2(ixx_fef_fefphase)<52);
fef_fefphase_ind_m2=(ind_ses_fef2(ixx_fef_fefphase)>51);

%% Fig 4a  FEF LFP IT Spike Correct versus Wrong
val_cr=spl_fef_it_cr(ixx_it_fefphase,:,:);val_wr=spl_fef_it_wr(ixx_it_fefphase,:,:);
bb1=0.2;bb2=0.3;
figure('name','Fig. 4a, Rezayat.E, et al')
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze((val_cr(:,1,:)));
var_h_wr=squeeze((val_wr(:,1,:)));
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),freqs2use,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),freqs2use,win_,c1(1),c1(2),c1(3))
xlabel('Frequency (Hz)')
ylabel('SPL (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
xlim([8 60])
ylim([bb1 bb2])
title('FEF LFP IT Spike')

export_file_to_excel([freqs2use',squeeze(nanmean(var_h_cr))',squeeze(nanmean(var_h_wr))'],'fig4a',13)

%% Fig 4a  FEF LFP IT Spike Correct versus Wrong
 ff=3;f=freqs2use;
f1=freq_bands(1,ff);  f2=freq_bands(2,ff); 
bb1=0.1; bb2=0.4;
    var_wr = [];var_wr = squeeze((val_wr(:,1,:)));
    var_cr = [];var_cr = squeeze((val_cr(:,1,:)));    
    x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));    
    y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2))); 
val_wr=x_h;val_cr=y_h;
ind_ses_m1=ixx_it_fefphase(it_fefphase_ind_m1);
ind_ses_m2=ixx_it_fefphase(it_fefphase_ind_m2);
figure('name','Fig. 4a inset, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr,val_wr);
n=length(x_h);
hold on
scatter(val_wr(it_fefphase_ind_m1),val_cr(it_fefphase_ind_m1),'b*')
scatter(val_wr(it_fefphase_ind_m2),val_cr(it_fefphase_ind_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[bb1 bb2];
ax.YTick=[bb1 bb2];
text(0.3,0.22,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(0.3,0.2,strcat('n=',num2str(n)),'fontsize',14);
xlim([bb1 bb2])
ylim([bb1 bb2])
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr,val_cr],'fig4a',1)

%% Fig 4b  IT LFP IT Spike Correct versus Wrong
val_cr=spl_it_it_cr(ixx_it_itphase,:,:);val_wr=spl_it_it_wr(ixx_it_itphase,:,:);
bb1=0.16;bb2=0.26;
figure('name','Fig. 4b, Rezayat.E, et al')
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze((val_cr(:,1,:)));
var_h_wr=squeeze((val_wr(:,1,:)));
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),freqs2use,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),freqs2use,win_,c1(1),c1(2),c1(3))
xlabel('Frequency (Hz)')
ylabel('SPL (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
xlim([8 60])
ylim([bb1 bb2])
title('IT LFP IT Spike')
export_file_to_excel([freqs2use',squeeze(nanmean(var_h_cr))',squeeze(nanmean(var_h_wr))'],'fig4b',13)

%% Fig 4b  IT LFP IT Spike Correct versus Wrong
bb1=0.1; bb2=0.4;
    var_wr = [];var_wr = squeeze((val_wr(:,1,:)));
    var_cr = [];var_cr = squeeze((val_cr(:,1,:)));    
    x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));    
    y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2))); 
val_wr=x_h;val_cr=y_h;
ind_ses_m1=ixx_it_itphase(it_itphase_ind_m1);
ind_ses_m2=ixx_it_itphase(it_itphase_ind_m2);
figure('name','Fig. 4b inset, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr,val_wr);
n=length(x_h);
hold on
scatter(val_wr(it_itphase_ind_m1),val_cr(it_itphase_ind_m1),'b*')
scatter(val_wr(it_itphase_ind_m2),val_cr(it_itphase_ind_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[bb1 bb2];
ax.YTick=[bb1 bb2];
text(0.3,0.22,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(0.3,0.2,strcat('n=',num2str(n)),'fontsize',14);
xlim([bb1 bb2])
ylim([bb1 bb2])
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr,val_cr],'fig4b',1)

%% Fig 4c  FEF LFP FEF Spike Correct versus Wrong
val_cr=spl_fef_fef_cr(ixx_fef_fefphase,:,:);val_wr=spl_fef_fef_wr(ixx_fef_fefphase,:,:);
bb1=0.15;bb2=0.25;
figure('name','Fig. 4c, Rezayat.E, et al')
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze((val_cr(:,1,:)));
var_h_wr=squeeze((val_wr(:,1,:)));
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),freqs2use,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),freqs2use,win_,c1(1),c1(2),c1(3))
xlabel('Frequency (Hz)')
ylabel('SPL (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
xlim([8 60])
ylim([bb1 bb2])
title('FEF LFP FEF Spike')
export_file_to_excel([freqs2use',squeeze(nanmean(var_h_cr))',squeeze(nanmean(var_h_wr))'],'fig4c',13)

%% Fig 4c  FEF LFP FEF Spike Correct versus Wrong
bb1=0.1; bb2=0.3;
    var_wr = [];var_wr = squeeze((val_wr(:,1,:)));
    var_cr = [];var_cr = squeeze((val_cr(:,1,:)));    
    x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));    
    y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2))); 
val_wr=x_h;val_cr=y_h;
ind_ses_m1=ixx_fef_fefphase(fef_fefphase_ind_m1);
ind_ses_m2=ixx_fef_fefphase(fef_fefphase_ind_m2);
figure('name','Fig. 4c inset, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr,val_wr);
n=length(x_h);
hold on
scatter(val_wr(fef_fefphase_ind_m1),val_cr(fef_fefphase_ind_m1),'b*')
scatter(val_wr(fef_fefphase_ind_m2),val_cr(fef_fefphase_ind_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[bb1 bb2];
ax.YTick=[bb1 bb2];
text(0.2,0.22,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(0.2,0.2,strcat('n=',num2str(n)),'fontsize',14);
xlim([bb1 bb2])
ylim([bb1 bb2])
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr,val_cr],'fig4c',1)

%% Fig 4d  IT LFP FEF Spike Correct versus Wrong
val_cr=spl_it_fef_cr(ixx_fef_itphase,:,:);val_wr=spl_it_fef_wr(ixx_fef_itphase,:,:);
bb1=0.21;bb2=0.31;
figure('name','Fig. 4d, Rezayat.E, et al')
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze((val_cr(:,1,:)));
var_h_wr=squeeze((val_wr(:,1,:)));
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),freqs2use,win_,c1(1),c1(2),c1(3))
c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),freqs2use,win_,c1(1),c1(2),c1(3))
xlabel('Frequency (Hz)')
ylabel('SPL (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
xlim([8 60])
ylim([bb1 bb2])
title('IT LFP FEF Spike')
export_file_to_excel([freqs2use',squeeze(nanmean(var_h_cr))',squeeze(nanmean(var_h_wr))'],'fig4d',13)

%% Fig 4d  IT LFP FEF Spike Correct versus Wrong
bb1=0.1; bb2=0.4;
    var_wr = [];var_wr = squeeze((val_wr(:,1,:)));
    var_cr = [];var_cr = squeeze((val_cr(:,1,:)));    
    x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));    
    y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2))); 
val_wr=x_h;val_cr=y_h;
ind_ses_m1=ixx_fef_itphase(fef_itphase_ind_m1);
ind_ses_m2=ixx_fef_itphase(fef_itphase_ind_m2);
figure('name','Fig. 4b inset, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr,val_wr);
n=length(x_h);
hold on
scatter(val_wr(fef_itphase_ind_m1),val_cr(fef_itphase_ind_m1),'b*')
scatter(val_wr(fef_itphase_ind_m2),val_cr(fef_itphase_ind_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[bb1 bb2];
ax.YTick=[bb1 bb2];
text(0.3,0.22,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(0.3,0.2,strcat('n=',num2str(n)),'fontsize',14);
xlim([bb1 bb2])
ylim([bb1 bb2])
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr,val_cr],'fig4d',1)

