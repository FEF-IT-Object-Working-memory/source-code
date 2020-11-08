function Fig1d_Power(resp)
%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 1d, Fig S3bc
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
freqs2use=[1:50 51:3:150 155:5:300];


t1_delay= 600; t2_delay= 1200;                                                % Delay period
t1_vis= 50; t2_vis= 350;                                                      % Visual period
t1_fix= -300; t2_fix= 50;                                                     % Fixation period
ind_b = [41 110];
st_ind=1; end_ind= 2100;

t = 1:2100;
t = decimate(t,5,'fir');
t_h = t-500;
similar_caxis = 1;bb1 = -3;bb2 = 3;bb1_ = -3; bb2_ = 3;

%% Load
gen_lfp_pairs(1).IT=resp(1).LFP_IT;
gen_lfp_pairs(1).FEF=resp(1).LFP_FEF;
gen_lfp_pairs(1).bhv{1}=resp(1).bhv;gen_lfp_pairs(1).condition{1}=resp(1).condition';
gen_lfp_pairs(1).preff_obj=resp(1).pref_obj_mu;
gen_lfp_pairs(1).npreff_obj=resp(1).npref_obj_mu;
gen_lfp_pairs(1).FEF_artifact=resp(1).FEF_artifact;
gen_lfp_pairs(1).IT_artifact=resp(1).IT_artifact;

%% Core calculation
disp('Wait! Calculating Power spectrum  **********')
cal=1;
if cal
    Power_it_cr=[];Power_it_wr=[];
    Power_fef_cr=[];Power_fef_wr=[];
    Power_it_cr_n=[];Power_it_wr_n=[];
    Power_fef_cr_n=[];Power_fef_wr_n=[];
    NS=size(gen_lfp_pairs,2);
    for ss=1:NS
        
        tic
        lfp=gen_lfp_pairs(ss);
        ix_h_cr=[];
        ix_h_wr=[];
        
        obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        for condi=1:6
            if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
            
            ix_h_cr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==1)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            ix_h_wr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==0)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            
        end
      
        %% Power
        
        s_fef_h = []; s_fef_h = lfp.FEF(:,st_ind:end_ind);
       
        s_it_h = []; s_it_h = lfp.IT(:,st_ind:end_ind);  
        
        p_it = []; p_fef = [];
        [p_it,p_fef,freqs2use,t]= wavelet_core_DS_pertrial(s_it_h,s_fef_h,freqs2use);
        P_it = []; P_it = abs(p_it).^2;
        P_fef = [];  P_fef = abs(p_fef).^2;
        
        sig_=[];sig_=P_it;
        sig_m_=squeeze(mean(sig_(:,:,ind_b(1):ind_b(2)),3));
        sig_m=mean(sig_m_);sig_m=repmat(sig_m,size(sig_,1),1,size(sig_,3));
        sig_std=std(sig_m_);sig_std=repmat(sig_std,size(sig_,1),1,size(sig_,3));
        P_it_n=(sig_-sig_m)./sig_std;
        
        sig_=[];sig_=P_fef;
        sig_m_=squeeze(mean(sig_(:,:,ind_b(1):ind_b(2)),3));
        sig_m=mean(sig_m_);sig_m=repmat(sig_m,size(sig_,1),1,size(sig_,3));
        sig_std=std(sig_m_);sig_std=repmat(sig_std,size(sig_,1),1,size(sig_,3));
        P_fef_n=(sig_-sig_m)./sig_std;
         
        for condi=1:2
            ix_cr_h_=(logical(ix_h_cr((condi-1)*3+1,:))|logical(ix_h_cr((condi-1)*3+2,:)))|logical(ix_h_cr((condi-1)*3+3,:));
            ix_wr_h_=(logical(ix_h_wr((condi-1)*3+1,:))|logical(ix_h_wr((condi-1)*3+2,:)))|logical(ix_h_wr((condi-1)*3+3,:));
            
            N_trials_cr = [];N_trials_cr =  sum(logical(ix_cr_h_));
            N_trials_wr = [];N_trials_wr =  sum(logical(ix_wr_h_));
            
            P_fef_cr = [];P_fef_cr = squeeze(P_fef(ix_cr_h_,:,:));
            Power_fef_cr(ss,condi,:,:)=squeeze(nanmean(P_fef_cr));
            P_fef_cr = [];P_fef_cr = squeeze(P_fef_n(ix_cr_h_,:,:));
            Power_fef_cr_n(ss,condi,:,:)=squeeze(nanmean(P_fef_cr));
            
            if N_trials_wr>1
                P_fef_wr = [];P_fef_wr = squeeze(P_fef(ix_wr_h_,:,:));
                Power_fef_wr(ss,condi,:,:)=squeeze(nanmean(P_fef_wr));
                P_fef_wr = [];P_fef_wr = squeeze(P_fef_n(ix_wr_h_,:,:));
                Power_fef_wr_n(ss,condi,:,:)=squeeze(nanmean(P_fef_wr));
            end
        end
        for condi=1:3
            ix_cr_h_=(logical(ix_h_cr(condi,:))|logical(ix_h_cr(condi+3,:)));
            ix_wr_h_=(logical(ix_h_wr(condi,:))|logical(ix_h_wr(condi+3,:)));
            
            N_trials_cr = [];N_trials_cr =  sum(logical(ix_cr_h_));
            N_trials_wr = [];N_trials_wr =  sum(logical(ix_wr_h_));
            P_it_cr = [];P_it_cr = squeeze(P_it(ix_cr_h_,:,:));
            Power_it_cr(ss,condi,:,:)=squeeze(nanmean(P_it_cr));
            P_it_cr = [];P_it_cr = squeeze(P_it_n(ix_cr_h_,:,:));
            Power_it_cr_n(ss,condi,:,:)=squeeze(nanmean(P_it_cr));
            
            if N_trials_wr>1
                P_it_wr = [];P_it_wr = squeeze(P_it(ix_wr_h_,:,:));
                Power_it_wr(ss,condi,:,:)=squeeze(nanmean(P_it_wr));
                P_it_wr = [];P_it_wr = squeeze(P_it_n(ix_wr_h_,:,:));
                Power_it_wr_n(ss,condi,:,:)=squeeze(nanmean(P_it_wr));
            end
            
            
        end
        
        
        toc
    end
%     save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Power\power_full.mat','t','freqs2use',...
%        'Power_it_cr','Power_it_wr','Power_fef_cr','Power_fef_wr',...
%        'Power_it_cr_n','Power_it_wr_n','Power_fef_cr_n','Power_fef_wr_n')
end

%% Load data
% load (strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Power\power_full.mat'))

% Baseline Correction
baseline_corrected=0;
if baseline_corrected
    plv_cr_=Power_it_cr;plv_wr_=Power_it_wr;
    val_b=nanmean(plv_cr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b,[],1);
    val_sd=repmat(val_sd,1,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_it_cr=(plv_cr_-val_b)./val_sd;
    
    val_b=nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b,[],1);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_it_wr=(plv_wr_-val_b)./val_sd;
    
    
    plv_cr_=Power_fef_cr;plv_wr_=Power_fef_wr;
    val_b=nanmean(plv_cr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b,[],1);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_fef_cr=(plv_cr_-val_b)./val_sd;
    
    val_b=nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b,[],1);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_fef_wr=(plv_wr_-val_b)./val_sd;
    
end

%% Smoothing For Hitmaps and Timecourse plots (Just for preview)
sigma_freq = 2;sigma_time = 15;
Power_fef_cr_smoothed=Power_fef_cr;Power_fef_wr_smoothed=Power_fef_wr;
for ss=1:size(Power_fef_cr_smoothed,1)
    for condi=1:2
        var_h=squeeze(Power_fef_cr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        Power_fef_cr_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(Power_fef_wr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        Power_fef_wr_smoothed(ss,condi,:,:)=var_h;      
        
      
    end
end

Power_it_cr_smoothed=Power_it_cr;Power_it_wr_smoothed=Power_it_wr;
for ss=1:size(Power_it_cr_smoothed,1)
    for condi=1:3
        var_h=squeeze(Power_it_cr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        Power_it_cr_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(Power_it_wr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        Power_it_wr_smoothed(ss,condi,:,:)=var_h;        
         
    end
end

%% Fig. 2d Power of FEF and IT for Cr and Wr trials
ind_ses=1;%session_selection(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
figure('name','Fig. 2d, Rezayat.E, et al')
for iteri=1:4
    ax= subplot(2,2,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (Power_fef_cr_smoothed);str = ' Cr';
        case 2
            P = (Power_fef_wr_smoothed);str = ' Wr';
        case 3
            P = (Power_it_cr_smoothed);str = ' Cr';
        case 4
            P = (Power_it_wr_smoothed);str = ' Wr';
    end
    if size(P,1)>1
        var_h=squeeze((P(ind_ses,1,:,:)));
        var_h=squeeze(nanmean((var_h),1));
    else
        var_h=squeeze((P(ind_ses,1,:,:)));
        var_h=squeeze(((var_h)));
    end
    h=pcolor(t_h/1000,freqs2use,var_h);
    h.EdgeColor = 'none';
    colormap(jet); clb= colorbar;
    ax.YDir='normal';
    ax.YScale='log';
    
    if similar_caxis;
        caxis([bb1 bb2]);%clb.Ticks=bb1:(bb2-bb1)/5:bb2;
        %         clb.TickLabels=num2cell(bb1:(bb2-bb1)/5:bb2);
        clb.Ticks=[-1 1];
    end
    if iteri==3
        clb.Label.String='Norm. PPL (a.u.)';
        
    end
    %     ax.YLim=[4 110];
    ax.XLim=[-0.300 1.500];
    line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',1)
    line([0.300 0.300],[ylim],'color','k','LineStyle',':','LineWidth',1)
    line([1.300 1.300],[ylim],'color','k','LineStyle',':','LineWidth',1)
    title(str)
    
    ax.XTick=[0 0.300 1.300];
    ax.YTick=[];
    if iteri==1|iteri==3
        ylabel('Frequency (Hz)')
       
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
     if iteri==3
     xlabel('Time from sample onset (sec.)')
     end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S4e Scatter Power Cr and Wr FEF
y_h= squeeze(nanmean(nanmean(Power_fef_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
x_h= squeeze(nanmean(nanmean(Power_fef_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));
ind_m1_=ind_ses_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_ses_m2((~isnan(x_h_)&~isnan(y_h_)));

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

%% State FEF Power  Cr vs. Wr 
try
disp('%%% State FEF Firing rate  Cr vs. Wr  %%%')

x_h=squeeze(nanmean(nanmean(Power_fef_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_fef_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(nanmean(Power_fef_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_fef_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(nanmean(Power_fef_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_fef_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('FEF In, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end

%% Fig. S4f Scatter Cr and Wr IT
y_h= squeeze(nanmean(nanmean(Power_it_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
x_h= squeeze(nanmean(nanmean(Power_it_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
x_h_=x_h;y_h_=y_h;
x_h = []; x_h = x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h = []; y_h = y_h_(~isnan(x_h_)&~isnan(y_h_));

ind_m1_=ind_ses_m1((~isnan(x_h_)&~isnan(y_h_)));
ind_m2_=ind_ses_m2((~isnan(x_h_)&~isnan(y_h_)));

p=signrank(x_h,y_h);
n=length(x_h);
 m_it_cr_wr_fr =  [nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h) ,n, p];

figure('name','Fig. S4f, Rezayat.E, et al')
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

%% State IT Power Cr vs. Wr 
try
disp('%%% State IT Firing rate  Cr vs. Wr  %%%')
x_h=squeeze(nanmean(nanmean(Power_it_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_it_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(nanmean(Power_it_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_it_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(nanmean(Power_it_wr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
y_h=squeeze(nanmean(nanmean(Power_it_cr(ind_ses,1,freqs2use>freq_bands(1,3)...
    &freqs2use<freq_bands(2,3),t1_delay<t_h&t_h<t2_delay),3),4));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('IT Pref, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end




end