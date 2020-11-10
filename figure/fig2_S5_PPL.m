%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 2, Fig S3

%% Opening 
clear all
clc
close all
disp ('??? ???? ?????? ?????? ')

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
t1_fix= -300; t2_fix= 50;                                                     % Fixation period

t = 1:2100;
t = decimate(t,5,'fir');
t_h = t-500;
similar_caxis = 1;bb1 = -1;bb2 = 1;bb1_ = -3; bb2_ = 3;

%% Load data

load([cd,'\Mat\PPL\PPL_shuff_1001_trial_matched.mat'])
ppl_cr_tri_matched = ppl_cr_n;ppl_wr_tri_matched = ppl_wr_n;

 load(strcat(cd,'\Mat\PPL\PPL_shuff_1001_2.mat'));
% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\PPL\PPL_full_shuff_1001.mat'));

% Shuffle Correction
ppl_cr_nsc=ppl_cr;ppl_wr_nsc=ppl_wr;
NS = size(ppl_cr_n,1);
shuffle_corrected = 1;
if shuffle_corrected
    ppl_cr = ppl_cr_n;ppl_wr = ppl_wr_n;
end

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

% Baseline Correction for not shuffle corretion data Fig. S3a
if baseline_corrected
    ppl_cr_h = ppl_cr_nsc;ppl_wr_h = ppl_wr_nsc;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr_nsc = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_wr_nsc = (ppl_wr_h-val_b)./val_sd;
end

% Baseline Correction for Trial matched data Fig. S3a
if baseline_corrected
    ppl_cr_h = ppl_cr_tri_matched;ppl_wr_h = ppl_wr_tri_matched;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr_tri_matched = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,1,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_wr_tri_matched = (ppl_wr_h-val_b)./val_sd;
end

%% Smoothing For Hitmaps and Timecourse plots (Just for preview)
sigma_freq = 2;sigma_time = 15;
ppl_cr_smoothed=ppl_cr;ppl_wr_smoothed=ppl_wr;
for ss=1:size(ppl_cr_smoothed,1)
    for condi=1:6
        var_h=squeeze(ppl_cr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_cr_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(ppl_wr(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_wr_smoothed(ss,condi,:,:)=var_h;
        
        
        var_h=squeeze(ppl_cr_nsc(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_cr_nsc_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(ppl_wr_nsc(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_wr_nsc_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(ppl_cr_tri_matched(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_cr_tri_match_smoothed(ss,condi,:,:)=var_h;
        
        var_h=squeeze(ppl_wr_tri_matched(ss,condi,:,:));
        var_h=imgaussfilt(var_h,[sigma_freq sigma_time]);
        ppl_wr_tri_match_smoothed(ss,condi,:,:)=var_h;
        
    end
end

%% Fig. 2a PPL between FEF and IT was different for Cr and Wr trials
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. 2a, Rezayat.E, et al')
for iteri=1:3
    ax= subplot(1,3,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (ppl_cr_smoothed);str = ' Cr';
        case 2
            P = (ppl_wr_smoothed);str = ' Wr';
        case 3
            P =( ppl_cr_smoothed)-(ppl_wr_smoothed);str = ' Cr - Wr';
    end
    var_h=squeeze((P(ind_ses,1,:,:)));
    var_h=squeeze(nanmean((var_h)));
%                   export_file_to_excel([var_h],['fig2a' str ] )
 x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'C' 'PJ' 'AFQ'};
xlswrite(file_name,x,'fig2a',strcat(col{iteri},'3'))
 
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. 2b Beta band PPL during sample and delay differed between Cr and Wr trials
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. 2b, Rezayat.E, et al')
f1=freq_bands(1,3);f2=freq_bands(2,3);
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze(nanmean(ppl_cr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3));
var_h_wr=squeeze(nanmean(ppl_wr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3));

c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),t_h/1000,1,c1(1),c1(2),c1(3))
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),t_h/1000,1,c1(1),c1(2),c1(3))
set(gca,'fontsize',14,'fontweight','bold')
              export_file_to_excel([t_h' nanmean(var_h_cr)' nanmean(var_h_wr)'],'fig2b',6 )
%               export_file_to_excel([],['fig2bwr'  ] )

xlabel('Time from sample onset (sec)');
ylabel('Norm. PPL (a.u)');

%
text(0.600,0.5,strcat(' Cr'),'fontsize',14,'fontweight','bold','color',Color_bhv(1,:));
text(0.600,-0.5,strcat(' Wr'),'fontsize',14,'fontweight','bold','color',Color_bhv(2,:));

ax.XLim = [-0.300 1.600]; ax.YLim=[-1 1];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[-1 1];
line([0 0],ylim,'Color','k');
line([0.300 0.300],ylim,'Color','k');
line([1.300 1.300],ylim,'Color','k');
line([0.600 1.200],[-0.08 -0.08],'Color','k','Linewidth',4);

%% Fig. 2c Comparison of beta band PPL during the delay period for correct versus wrong trials
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. 2c, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);

export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'fig2c',1)

%% Fig. 2c Comparison of beta band PPL during the delay period for correct versus wrong trials
bin_=0.5;
hist_var=[];nbin1=-2:bin_:2.5;nbin2=-2.5:bin_:2;

figure('name','Fig. 2c, Rezayat.E, et al')

subplot(211)
hold on
[counts,centers]=hist(val_wr(ind_ses,1),nbin1);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)
xlim([-2 3])
m=nanmedian(val_wr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_wr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_wr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

subplot(212)
hold on
[counts,centers]=hist(val_cr(ind_ses),nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
xlim([-3 2])
line([0 0],ylim)
m=nanmedian(val_cr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_cr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_cr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

%% Fig. 2d Comparison of object selectivity of beta band PPL during the delay period for correct versus wrong trials
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))-...
    squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))-...
    squeeze(nanmean(nanmean(ppl_wr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. 2d, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('PPL object sel. (a.u.) [Cr]')
xlabel('PPL object sel. (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'fig2d',1)

%% Fig. 2d Comparison of object selectivity of beta band PPL during the delay period for correct versus wrong
bin_=0.5;
hist_var=[];nbin1=-2:bin_:2.5;nbin2=-2.5:bin_:2;

figure('name','Fig. 2d, Rezayat.E, et al')

subplot(211)
hold on
[counts,centers]=hist(val_wr(ind_ses,1),nbin1);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)
xlim([-2 3])
m=nanmedian(val_wr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_wr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_wr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

subplot(212)
hold on
[counts,centers]=hist(val_cr(ind_ses),nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
xlim([-3 2])
line([0 0],ylim)
m=nanmedian(val_cr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_cr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_cr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

%% Fig. 2e Comparison of location selectivity of beta band PPL during the delay period for correct versus wrong trials
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))-...
    squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))-...
    squeeze(nanmean(nanmean(ppl_wr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. 2e, Rezayat.E, et al')
ax=subplot(1,1,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('PPL location sel. (a.u.) [Cr]')
xlabel('PPL location sel. (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
%               export_file_to_excel([val_wr val_cr],['fig2e'  ] )
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'fig2e',1)

%% Fig. 2e Comparison of location selectivity of beta band PPL during the delay period for correct versus wrong
bin_=0.5;
hist_var=[];nbin1=-2:bin_:2.5;nbin2=-2.5:bin_:2;

figure('name','Fig. 2e, Rezayat.E, et al')

subplot(211)
hold on
[counts,centers]=hist(val_wr(ind_ses,1),nbin1);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)
xlim([-2 3])
m=nanmedian(val_wr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_wr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_wr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

subplot(212)
hold on
[counts,centers]=hist(val_cr(ind_ses),nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
xlim([-3 2])
line([0 0],ylim)
m=nanmedian(val_cr(ind_ses));
line([m m],ylim,'linewidth',2,'color','k')

m1=nanmedian(val_cr(ind_ses_m1));
line([m1 m1],ylim,'linewidth',2,'color','b')

m2=nanmedian(val_cr(ind_ses_m2));
line([m2 m2],ylim,'linewidth',2,'color','g')

%% Fig. S5a Beta band PPL was greater for correct trials without the shuffle correction
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5a, Rezayat.E, et al')
for iteri=1:2
    ax= subplot(1,2,iteri);
    str=[];P = [];
     str=[];P = [];
    switch iteri
        case 1
         P =( ppl_cr_tri_match_smoothed)-(ppl_wr_tri_match_smoothed);str = ' Cr - Wr [Trial matched]';
        case 2
         P =( ppl_cr_nsc_smoothed)-(ppl_wr_nsc_smoothed);str = ' Cr - Wr [Not Shuffled]';
       
    end
    
    var_h=squeeze((P(ind_ses,1,:,:)));
    var_h=squeeze(nanmean((var_h)));
        x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'C' 'PJ'};
xlswrite(file_name,x,'figS5a',strcat(col{iteri},'3'))
           
%     export_file_to_excel([var_h ],['fig5Sa' str ] )

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
    if iteri==2
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end
%% Fig. S5a Beta band PPL was greater for correct trials without the shuffle correction
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5a, Rezayat.E, et al')
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr_tri_matched(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_wr_tri_matched(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
ax=subplot(2,1,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
title('Trial matched')
%               export_file_to_excel([val_wr val_cr],['figS5amatch'  ] )
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5a_match',1)

val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr_nsc(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_wr_nsc(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
ax=subplot(2,1,2);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Cr]')
xlabel('Norm. PPL (a.u.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
title('Not shuffled')
%               export_file_to_excel([val_wr val_cr],['figS5ashuff'  ] )
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5a_shuff',1)

%% Fig. S5b The trend toward higher beta band PPL on correct trials was present in both M1 and M2
figure('name','Fig. S5b, Rezayat.E, et al')
for iteri=1:2
    ax= subplot(1,2,iteri);
    str=[];
    ix_h = []; if ismember(iteri,[1 3]);ix_h = ind_ses_m1;str='M1'; else ix_h = ind_ses_m2; str='M2'; end
    P = []; P = ppl_cr_smoothed-ppl_wr_smoothed;str=[str ''];
    var_h=squeeze((P(ix_h,1,:,:)));
    var_h=squeeze(nanmean((var_h)));
    h=pcolor(t_h/1000,freqs2use,var_h);
%                           export_file_to_excel([var_h ],['figS5b' str ] )
 x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'C' 'PJ'};
xlswrite(file_name,x,'figS5b',strcat(col{iteri},'3'))
 
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
    if iteri==2
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
    end
    text(0.5,5,strcat('n=',num2str(length(ix_h))),'fontsize',14)
    
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S5c Inter-areal beta PPL encoded object identity during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5c, Rezayat.E, et al')
for iteri=1:3
    ax= subplot(1,3,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (ppl_cr_smoothed(ind_ses,1,:,:));str = ' Pref';
        case 2
            P = (ppl_cr_smoothed(ind_ses,3,:,:));str = ' NPref';
        case 3
            P =( ppl_cr_smoothed(ind_ses,1,:,:))-(ppl_cr_smoothed(ind_ses,3,:,:));str = ' Pref - NPref';
    end
    var_h=squeeze((P));
    var_h=squeeze(nanmean((var_h)));
%                       export_file_to_excel([var_h ],['figS5c' str ] )
x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'D' 'PJ' 'AFQ'};
xlswrite(file_name,x,'figS5c',strcat(col{iteri},'3'))
 
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S5c Inter-areal beta PPL encoded object identity during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5c, Rezayat.E, et al')
f1=freq_bands(1,3);f2=freq_bands(2,3);
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze(nanmean(ppl_cr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3))-...
    squeeze(nanmean(ppl_cr_smoothed(ind_ses,3,freqs2use>f1&freqs2use<f2,:),3));
var_h_wr=squeeze(nanmean(ppl_wr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3))-...
    squeeze(nanmean(ppl_wr_smoothed(ind_ses,3,freqs2use>f1&freqs2use<f2,:),3));

c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),t_h/1000,1,c1(1),c1(2),c1(3))
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),t_h/1000,1,c1(1),c1(2),c1(3))
set(gca,'fontsize',14,'fontweight','bold')

xlabel('Time from sample onset (sec)');
ylabel('PPL object sel. (a.u)');

%
text(0.600,0.5,strcat(' Cr'),'fontsize',14,'fontweight','bold','color',Color_bhv(1,:));
text(0.600,-0.5,strcat(' Wr'),'fontsize',14,'fontweight','bold','color',Color_bhv(2,:));

ax.XLim = [-0.300 1.600]; ax.YLim=[-1 1];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[-1 1];
line([0 0],ylim,'Color','k');
line([0.300 0.300],ylim,'Color','k');
line([1.300 1.300],ylim,'Color','k');
line([0.600 1.200],[-0.08 -0.08],'Color','k','Linewidth',4);
x=[t_h' nanmean(var_h_cr)' nanmean(var_h_wr)'];
a=[]; a={'time','Obj PPL Cr','Obj PLL Wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,'figS5c','B85')

%% Fig. S5c Inter-areal beta PPL encoded object identity during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
f1=freq_bands(1,3);f2=freq_bands(2,3);
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. S5c, Rezayat.E, et al')
ax=subplot(1,2,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Pref]')
xlabel('Norm. PPL (a.u.) [NPref]')
title('\beta')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5cscbeta',1)

f1=freq_bands(1,4);f2=freq_bands(2,4);

val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
ax=subplot(1,2,2);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [Pref]')
xlabel('Norm. PPL (a.u.) [NPref]')
title('L \gamma')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
 export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5cscgamma',1)

%% Fig. S5d Inter-areal beta PPL encoded the location of the sample during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5d, Rezayat.E, et al')
for iteri=1:3
    ax= subplot(1,3,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (ppl_cr_smoothed(ind_ses,1,:,:));str = ' In';
        case 2
            P = (ppl_cr_smoothed(ind_ses,4,:,:));str = ' Out';
        case 3
            P =( ppl_cr_smoothed(ind_ses,1,:,:))-(ppl_cr_smoothed(ind_ses,4,:,:));str = ' In - Out';
    end
    var_h=squeeze((P));
    var_h=squeeze(nanmean((var_h)));
%                       export_file_to_excel([var_h ],['figS5d' str ] )
x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'D' 'PJ' 'AFQ'};
xlswrite(file_name,x,'figS5d',strcat(col{iteri},'3'))
 
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S5d Inter-areal beta PPL encoded the location of the sample during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S5c, Rezayat.E, et al')
f1=freq_bands(1,3);f2=freq_bands(2,3);
ax= subplot(1,1,1);
hold on
var_h_cr=squeeze(nanmean(ppl_cr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3))-...
    squeeze(nanmean(ppl_cr_smoothed(ind_ses,4,freqs2use>f1&freqs2use<f2,:),3));
var_h_wr=squeeze(nanmean(ppl_wr_smoothed(ind_ses,1,freqs2use>f1&freqs2use<f2,:),3))-...
    squeeze(nanmean(ppl_wr_smoothed(ind_ses,4,freqs2use>f1&freqs2use<f2,:),3));

c1 = Color_bhv(1,:);
niceplot(squeeze((var_h_cr)),t_h/1000,1,c1(1),c1(2),c1(3))
c1 = Color_bhv(2,:);
niceplot(squeeze((var_h_wr)),t_h/1000,1,c1(1),c1(2),c1(3))
set(gca,'fontsize',14,'fontweight','bold')

xlabel('Time from sample onset (sec)');
ylabel('PPL location sel. (a.u)');

%
text(0.600,0.5,strcat(' Cr'),'fontsize',14,'fontweight','bold','color',Color_bhv(1,:));
text(0.600,-0.5,strcat(' Wr'),'fontsize',14,'fontweight','bold','color',Color_bhv(2,:));

ax.XLim = [-0.300 1.600]; ax.YLim=[-1 1];
set(gca,'fontsize',14,'fontweight','bold');
ax.XTick=[0 0.300 1.300];
ax.YTick=[-1 1];
line([0 0],ylim,'Color','k');
line([0.300 0.300],ylim,'Color','k');
line([1.300 1.300],ylim,'Color','k');
line([0.600 1.200],[-0.08 -0.08],'Color','k','Linewidth',4);

x=[t_h' nanmean(var_h_cr)' nanmean(var_h_wr)'];
a=[]; a={'time','Loc PPL Cr','Loc PLL Wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,'figS5d','B85')

%% Fig. S5d Inter-areal beta PPL encoded the location of the sample during the delay period
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
f1=freq_bands(1,3);f2=freq_bands(2,3);
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. S5c, Rezayat.E, et al')
ax=subplot(1,2,1);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [In]')
xlabel('Norm. PPL (a.u.) [Out]')
title('\beta')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5dscbeta',1)

f1=freq_bands(1,4);f2=freq_bands(2,4);

val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
ax=subplot(1,2,2);
p=signrank(val_cr(ind_ses),val_wr(ind_ses));
n=length(ind_ses);
hold on
scatter(val_wr(ind_ses_m1),val_cr(ind_ses_m1),'b*')
scatter(val_wr(ind_ses_m2),val_cr(ind_ses_m2),'gs')
plot([-2:0.1:5],[-2:0.1:5],'k')
axis([-2 3 -2 3])
ylabel('Norm. PPL (a.u.) [In]')
xlabel('Norm. PPL (a.u.) [Out]')
title('L \gamma')
set(gca,'fontsize',14,'fontweight','bold')
ax.XTick=[-2 3];
ax.YTick=[-2 3];
text(2,-1,strcat('p=',num2str(round(p,3))),'fontsize',14);
text(2,-1.5,strcat('n=',num2str(n)),'fontsize',14);
 export_file_to_excel([[ones(length(ind_ses_m1),1); 2*ones(length(ind_ses_m2),1)],val_wr(ind_ses),val_cr(ind_ses)],'figS5dscgamma',1)

%% State for beta-band PPL during the delay period compared to baseline For Correct trials
disp('%%% State for beta-band PPL during the delay period compared to baseline For Correct trials %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_fix<t_h&t_h<t2_fix),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
corretc_increase=sum((y_h-x_h)>0);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; ',...
    num2str(corretc_increase),'sessions show increase in Cr trials'));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

%% State for beta-band PPL during the delay period compared to baseline For Wrong trials
disp('%%% State for beta-band PPL during the delay period compared to baseline For Wrong trials %%%')
ind_ses=session_selection_PPL(3);                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_fix<t_h&t_h<t2_fix),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
wrong_decrease=sum((y_h-x_h)<0);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; ',...
    num2str(wrong_decrease),'sessions show decrease in Wr trials'));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

%% State for beta-band PPL during the delay period for Correct trials
disp('%%% State for beta-band PPL during the delay period for Correct trials %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_fix<t_h&t_h<t2_fix),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Delay vs. Baseline, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for object selectivity of beta-band PPL during the delay period
disp('%%% State for object selectivity of beta-band PPL during the Delay period %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Pref vs. NPref, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Pref vs. NPref, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Pref vs. NPref, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for location selectivity of beta-band PPL during the delay period
disp('%%% State for location selectivity of beta-band PPL during the Delay period %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, In vs. Out, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, In vs. Out, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Behavioral correlation of beta-band PPL during the delay period
disp('%%% State for Behavioral correlation of beta-band PPL during the Delay period %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Behavioral correlation of object selectivity during the Delay period
disp('%%% State for Behavioral correlation of object selectivity during the Delay period %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))...
    -squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))...
    -squeeze(nanmean(nanmean(ppl_wr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Object sel. PPL Delay, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Object sel.PPL Delay, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Object sel.PPL Delay, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Behavioral correlation of location selectivity during the Delay period
disp('%%% State for Behavioral correlation of location selectivity during the Delay period %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))...
    -squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3))...
    -squeeze(nanmean(nanmean(ppl_wr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Location sel. PPL Delay, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Location sel. PPL Delay, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Location sel. PPL Delay, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Trial matched]
disp('%%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Trial matched] %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr_tri_matched(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr_tri_matched(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Not Shuffle Corrected]
disp('%%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Not Shuffle Corrected] %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr_nsc(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr_nsc(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Delay, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

%% State for beta-band PPL during the visual period for Correct trials
disp('%%% State for beta-band PPL during the visual period for Correct trials %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_fix<t_h&t_h<t2_fix),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Visual vs. Baseline, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Visual vs. Baseline, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL, Visual vs. Baseline, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for object selectivity of beta-band PPL during the visual period
disp('%%% State for object selectivity of beta-band PPL during the visual period %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Pref vs. NPref, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Pref vs. NPref, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Pref vs. NPref, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for location selectivity of beta-band PPL during the visual period
disp('%%% State for location selectivity of beta-band PPL during the visual period %%%')
ind_ses=session_selection_PPL(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, In vs. Out, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, In vs. Out, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
%% State for Behavioral correlation of beta-band PPL during the visual period
disp('%%% State for Behavioral correlation of beta-band PPL during the visual period %%%')
ind_ses=session_selection_PPL(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

f1=freq_bands(1,3);f2=freq_bands(2,3);
val_y_h=[];val_y_h=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));
val_x_h=[];val_x_h=squeeze(nanmean(nanmean(ppl_wr(:,1,freqs2use>f1&freqs2use<f2,t1_vis<t_h&t_h<t2_vis),4),3));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses);y_h=val_y_h(ind_ses);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),'; '));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m1);y_h=val_y_h(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=[];y_h=[];
x_h=val_x_h(ind_ses_m2);y_h=val_y_h(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('PPL Visual, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));


