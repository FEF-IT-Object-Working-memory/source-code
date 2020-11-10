%% In the name of Allah
%Power spectrum wavelet versions 2
clear all
clc
close all

%% Settings
obj_label={'Pref' 'Intermediate' 'NPref'};
loc_label={'IN' 'OUT'};
N_shuff=100;
func=@ (x) (x);
g_f=2;g_t=15;
ind_b=[41 110];

freq_bands=[3.9,7.9,21.9,34.9,49.9,150;8.1,15.1,28.1,50.1,130.1,300];

bands_h_l={ ' \theta' ' \alpha' ' \beta' '  L \gamma' ' H \gamma' ' 150< \gamma'};
bands_h_name={ 'theta' 'alpha' 'beta' 'Lgamma'  'Hgamma' '150gamma'};

t1= 600; t2= 1200;

%% Load
load ([ cd ,'\Mat\Power\power_full.mat'])
NS=size(Power_it_cr,1); 
t_h=t-500;
f=freqs2use;
similar_caxis=1;bb1=-1.5;bb2=1.5; bb1_=-0.2;bb2_=0.2;
 bb1_=bb1;bb2_=bb2;

 %%
 baseline_corrected=1;
if baseline_corrected
    plv_cr_=Power_it_cr;plv_wr_=Power_it_wr;
    val_b=nanmean(plv_cr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_it_cr=(plv_cr_-val_b)./val_sd;
    
    val_b=nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_it_wr=(plv_wr_-val_b)./val_sd;
    
    
    plv_cr_=Power_fef_cr;plv_wr_=Power_fef_wr;
    val_b=nanmean(plv_cr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_fef_cr=(plv_cr_-val_b)./val_sd;
    
    val_b=nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd=nanstd(val_b);
    val_sd=repmat(val_sd,NS,1,1);
    val_sd=repmat(val_sd,1,1,1,420);
    val_b=repmat(val_b,1,1,1,420);
    Power_fef_wr=(plv_wr_-val_b)./val_sd;
    
end

%% Session selection for comparing corect and Wrong

     load multiunit_preferance_pval.mat
     load('prefereance_power.mat')
     p=0.05;
     ind_nonslective=find((p_rate>p)&(p_baseline>p));


ind_ses_fef=session_select_power_after_review(8);
ind_m1_fef=(ind_ses_fef<52);
ind_m2_fef=(ind_ses_fef>51);

ind_ses_it=session_select_power_after_review(10);
ind_ses_it=setdiff(ind_ses_it,ind_nonslective);
ind_m1_it=(ind_ses_it<52);
ind_m2_it=(ind_ses_it>51);

ind_ses_m1_fef=setdiff(ind_ses_fef,[52:NS]);
ind_ses_m1_it=setdiff(ind_ses_it,[52:NS]);
ind_ses_m2_fef=setdiff(ind_ses_fef,[1:51]);
ind_ses_m2_it=setdiff(ind_ses_it,[1:51]);

%% Fig 1d Preview Correct versus Wrong FEF
figure('name','CR & WR IN Pref')
for iteri=1:2
    ax= subplot(1,2,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = func(Power_fef_cr(ind_ses_fef,1,:,:));str = 'Correct';
        case 2
            P = func(Power_fef_wr(ind_ses_fef,1,:,:));str = 'Wrong';
        case 3
            P = func(Power_fef_cr(ind_ses_fef,1,:,:))-func(Power_fef_wr(ind_ses_fef,1,:,:));str = 'Diff';
    end

    var_h=squeeze((P));
    var_h=squeeze(nanmean((var_h)));
    var_h=imgaussfilt(var_h,[g_f g_t]);
   x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'C' 'PJ'};
xlswrite(file_name,x,'fig1d',strcat(col{iteri},'3'))
    

    h=pcolor(t_h,freqs2use,var_h);
    h.EdgeColor = 'none';
    colormap(jet); clb= colorbar;
    ax.YDir='normal';
    ax.YScale='log';
    ax.YTick=[4  15 35  100 ];
    ax.XTick=[0 300 1300];
    
    if similar_caxis;
        if iteri==3|iteri==6
            caxis([bb1_ bb2_]);clb.Ticks=bb1_:(bb2_-bb1_)/5:bb2_;
            clb.TickLabels=num2cell(bb1_:(bb2_-bb1_)/5:bb2_);
            clb.TickLabels=num2cell(bb1_:(bb2_-bb1_)/5:bb2_);
            clb.Label.String='Baseline normalize power (a.u)';
            line([600 1200],[4 4],'color','k','LineWidth',4)
        else
            caxis([bb1 bb2]);clb.Ticks=bb1:(bb2-bb1)/5:bb2;
            clb.TickLabels=num2cell(bb1:(bb2-bb1)/5:bb2);
            clb.TickLabels=num2cell(bb1:(bb2-bb1)/5:bb2);
        end
    end
    
%     ax.YLim=[4 150];ax.XLim=[-300 1500];
    line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
    line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
    line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
    title(str)
    xlabel('Time from sample onset (ms)')
    if iteri==1| iteri==4
        ylabel('Frequency (Hz)')
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig S4e Scatter Cr versus Wr IN  FEF
% ind_ses_fef=setdiff(ind_ses_fef,[52:62]);
condi=1;
figure('name','FEF Cr&Wr IN ')
 ff=3;
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    ax=subplot(1,1,1);
    hold on
    
    var_wr = [];var_wr = squeeze((Power_fef_wr(ind_ses_fef,condi,:,:)));
    var_cr = [];var_cr = squeeze((Power_fef_cr(ind_ses_fef,condi,:,:)));
    
    x_h=[]; x_h=func(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    y_h=[];y_h=func(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
        
    
    p=signrank(x_h,y_h);    n=length(x_h);
    eval(strcat('m_fef_cr_wr', bands_h_name{ff},'= ([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n) ,n, p])'))
    
    disp (strcat('Power FEF  Cr vs. Wr, All;, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

    
    hold on
%     y_h=func(squeeze(nanmean(nanmean(Power_fef_cr(ind_ses_fef_si,1,f>f1&f<f2,t1<t_h&t_h<t2),4),3)));
%     x_h=func(squeeze(nanmean(nanmean(Power_fef_wr(ind_ses_fef_si,1,f>f1&f<f2,t1<t_h&t_h<t2),4),3)));
%     scatter(x_h,y_h,'k','filled')
%     y_h=func(squeeze(nanmean(nanmean(Power_fef_cr(ind_ses_fef_nsi,1,f>f1&f<f2,t1<t_h&t_h<t2),4),3)));
%     x_h=func(squeeze(nanmean(nanmean(Power_fef_wr(ind_ses_fef_nsi,1,f>f1&f<f2,t1<t_h&t_h<t2),4),3)));
%     scatter(x_h,y_h,'k')
    
    % scatter(x_h,y_h,'k')
     scatter(x_h(ind_m1_fef),y_h(ind_m1_fef),'b*')
scatter(x_h(ind_m2_fef),y_h(ind_m2_fef),'gs')

     
    if p>0.00099
        text( 0.05,-0.05,strcat('p = ',num2str(p,3)),'fontsize',14,'fontweight','bold')
        
    else
        text( 0.05,-0.1,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
    end
    text( 0.05,-0.3,strcat('n = ',num2str(n,2)),'fontsize',14,'fontweight','bold')
%           export_file_to_excel([x_h y_h],['figS4e' ] )
export_file_to_excel([[ones(sum(ind_m1_fef),1); 2*ones(sum(ind_m2_fef),1)],x_h,y_h],'figS4e',1)

    plot([-2:0.1:10],[-2:0.1:10],'k:')
    set(gca,'fontsize',14,'fontweight','bold')
    if ff>2;    xlabel('Wr');end
    if mod(ff,2);    ylabel('Cr');end
    title(strcat(bands_h_l{ff}))
    
    axis([-1 1 -1 1])
  
%% Fig 1d Preview Correct versus Wrong IT
figure('name','CR & WR IN Pref')
for iteri=1:2
    ax= subplot(1,2,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = func(Power_it_cr(ind_ses_it,1,:,:));str = 'Correct';
        case 2
            P = func(Power_it_wr(ind_ses_it,1,:,:));str = 'Wrong';
        case 3
            P =func( Power_it_cr(ind_ses_it,1,:,:))-func(Power_it_wr(ind_ses_it,1,:,:));str = 'Diff';
    end
    var_h=squeeze((P));
    var_h=squeeze(nanmean((var_h)));
    var_h=imgaussfilt(var_h,[g_f g_t]);
%                   export_file_to_excel([var_h],['fig1d_it' str ] )
   x=[t_h;var_h];
   x=[[0 freqs2use]',x];
    file_name=strcat('C:\Users\Ehsan\Desktop\PDF\filedata_final.xlsx');
col={'C' 'PJ'};
xlswrite(file_name,x,'fig1d',strcat(col{iteri},'122'))
 
    h=pcolor(t_h,freqs2use,var_h);
    h.EdgeColor = 'none';
    colormap(jet); clb= colorbar;
    ax.YDir='normal';
    ax.YScale='log';
    ax.YTick=[4  15 35  100 ];
    ax.XTick=[0 300 1300];
    
    if similar_caxis;
        if iteri==3|iteri==6
            caxis([bb1_ bb2_]);clb.Ticks=bb1_:(bb2_-bb1_)/5:bb2_;
            clb.TickLabels=num2cell(bb1_:(bb2_-bb1_)/5:bb2_);
            clb.TickLabels=num2cell(bb1_:(bb2_-bb1_)/5:bb2_);
            clb.Label.String='Baseline normalize power (a.u)';
            line([600 1200],[4 4],'color','k','LineWidth',4)
        else
            caxis([bb1 bb2]);clb.Ticks=bb1:(bb2-bb1)/5:bb2;
            clb.TickLabels=num2cell(bb1:(bb2-bb1)/5:bb2);
            clb.TickLabels=num2cell(bb1:(bb2-bb1)/5:bb2);
        end
    end
    
%     ax.YLim=[4 150];ax.XLim=[-300 1500];
    line([0 0],[ylim],'color','k','LineStyle',':','LineWidth',2)
    line([300 300],[ylim],'color','k','LineStyle',':','LineWidth',2)
    line([1300 1300],[ylim],'color','k','LineStyle',':','LineWidth',2)
    
    
    if condi<4;loc=1;obj=condi;else loc=2;obj=condi-3;end
    title(str)
    
    xlabel('Time from sample onset (ms)')
    
    if iteri==1| iteri==4
        ylabel('Frequency (Hz)')
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig S4f Scatter Cr versus Wr Pref  IT
condi=1;
figure('name','IT Cr &Wr Preff')
ff=3;
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    ax=subplot(1,1,1);
    hold on
    var_wr = [];var_wr = squeeze((Power_it_wr(ind_ses_it,condi,:,:)));
    var_cr = [];var_cr = squeeze((Power_it_cr(ind_ses_it,condi,:,:)));
    x_h=func(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    y_h=func(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    
    %
    p=signrank(x_h,y_h);  n=length(x_h);
    
    eval(strcat('m_it_cr_wr', bands_h_name{ff},'= ([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p])'))
     disp (strcat('Power IT  Cr vs. Wr, All;, Effect size',' =',num2str( ...
       ((nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h)))/nanstd(y_h-x_h))))

      scatter(x_h(ind_m1_it),y_h(ind_m1_it),'b*')
scatter(x_h(ind_m2_it),y_h(ind_m2_it),'gs')
 
%     scatter(x_h,y_h,'k')
    if p>0.00099
        text( 0.05,-0.05,strcat('p = ',num2str(p,3)),'fontsize',14,'fontweight','bold')
        
    else
        text( 0.05,-0.1,strcat('p < 0.001'),'fontsize',14,'fontweight','bold')
    end
    text( 0.05,-0.3,strcat('n = ',num2str(length(x_h),2)),'fontsize',14,'fontweight','bold')
%       export_file_to_excel([x_h y_h],['figS4f' ] )
export_file_to_excel([[ones(sum(ind_m1_it),1); 2*ones(sum(ind_m2_it),1)],x_h,y_h],'figS4f',1)

    plot([-20:0.1:50],[-20:0.1:50],'k:')
    set(gca,'fontsize',14,'fontweight','bold')
    if ff>2;    xlabel('Wr');end
    if mod(ff,2);    ylabel('Cr');end
    title(strcat(bands_h_l{ff}))
    axis([-1 1 -1 1])
    
%%
Color_loc=[1,0,0;0,0,1];
Color_obj=[1,0.5,0;0.5,0.5,0.5;0,1,0.5];
Color_bhv=[1,1,0;0,1,1];
win_=10;

%% Late delay response stat test
%%% FEF  Cr Wr
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_fef_cr(ind_ses_fef,1,:,:));
    var_cr4=squeeze(Power_fef_wr(ind_ses_fef,1,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
  
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_fef_del_crwr', bands_h_name{ff},'= ([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n) ,n, p1])'))

  ind_h_m=(ind_ses_fef<52);
     p1=signrank(x_h(ind_h_m),y_h(ind_h_m))
    n=length(x_h(ind_h_m));    
    eval(strcat('m_fef_del_crwr_m1', bands_h_name{ff},'= ([nanmean(y_h(ind_h_m)) , nanmean(x_h(ind_h_m)), nanmean(y_h(ind_h_m))-nanmean(x_h(ind_h_m)),nanstd(y_h(ind_h_m)-x_h(ind_h_m))/sqrt(n) ,n, p1])'))

     ind_h_m=(ind_ses_fef>51);
     p1=signrank(x_h(ind_h_m),y_h(ind_h_m))
    n=length(x_h(ind_h_m));    
    eval(strcat('m_fef_del_crwr_m2', bands_h_name{ff},'= ([nanmean(y_h(ind_h_m)) , nanmean(x_h(ind_h_m)), nanmean(y_h(ind_h_m))-nanmean(x_h(ind_h_m)),nanstd(y_h(ind_h_m)-x_h(ind_h_m))/sqrt(n) ,n, p1])'))

end

%%% IT Cr Wr
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_it_cr(ind_ses_it,1,:,:));
    var_cr4=squeeze(Power_it_wr(ind_ses_it,1,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
   
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_it_del_crwr', bands_h_name{ff},'= ([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n) ,n, p1])'))

    ind_h_m=(ind_ses_it<52);
     p1=signrank(x_h(ind_h_m),y_h(ind_h_m))
    n=length(x_h(ind_h_m));    
    eval(strcat('m_it_del_crwr_m1', bands_h_name{ff},'= ([nanmean(y_h(ind_h_m)) , nanmean(x_h(ind_h_m)), nanmean(y_h(ind_h_m))-nanmean(x_h(ind_h_m)),nanstd(y_h(ind_h_m)-x_h(ind_h_m))/sqrt(n) ,n, p1])'))

     ind_h_m=(ind_ses_it>51);
     p1=signrank(x_h(ind_h_m),y_h(ind_h_m))
    n=length(x_h(ind_h_m));    
    eval(strcat('m_it_del_crwr_m2', bands_h_name{ff},'= ([nanmean(y_h(ind_h_m)) , nanmean(x_h(ind_h_m)), nanmean(y_h(ind_h_m))-nanmean(x_h(ind_h_m)),nanstd(y_h(ind_h_m)-x_h(ind_h_m))/sqrt(n) ,n, p1])'))

end

% %%% Gamma selective FEF Cr Wr
% for ff=1:size(freq_bands,2)
%     
%     f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
%     
%     var_cr1=squeeze(Power_fef_cr(ind_ses_fef_si,1,:,:));
%     var_cr4=squeeze(Power_fef_wr(ind_ses_fef_si,1,:,:));
%     
%     x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
%    
%     y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
%     
%     p1=signrank(x_h,y_h)
%     n=length(x_h);    
%     eval(strcat('m_fef_del_crwr_selec', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h) ,n, p1],2)'))
% 
%   
% end

%% Session selection for corect 
ind_ses_fef=session_select_power_after_review(1);
ind_m1_fef=(ind_ses_fef<52);
ind_m2_fef=(ind_ses_fef>51);
ind_ses_it=session_select_power_after_review(7);
 ind_ses_it=setdiff(ind_ses_it,ind_nonslective);
ind_m1_it=(ind_ses_it<52);
ind_m2_it=(ind_ses_it>51);
ind_ses_m1_fef=setdiff(ind_ses_fef,[52:NS]);
ind_ses_m1_it=setdiff(ind_ses_it,[52:NS]);
ind_ses_m2_fef=setdiff(ind_ses_fef,[1:51]);
ind_ses_m2_it=setdiff(ind_ses_it,[1:51]);

% xx=(pval_dsi_fef<0.05);
%  ind_ses_fef_si=find(sum(xx(:,freqs2use>35),2)>5);
% % ind_ses_fef_si=find(sum(xx(:,freqs2use>18&freqs2use<29),2)>5);
% 
% ind_ses_fef_nsi=setdiff(ind_ses_fef,ind_ses_fef_si);
% xx=(pval_dsi_it<0.05);
% ind_ses_it_si=find(sum(xx(:,freqs2use>8&freqs2use<15),2)>2);
% ind_ses_it_nsi=setdiff(ind_ses_it,ind_ses_it_si);

%% Visual response stat test
%%% FEF IN & OUT Correct
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_fef_cr(ind_ses_fef,1,:,:));
    var_cr4=squeeze(Power_fef_cr(ind_ses_fef,2,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,50<t_h&t_h<350),3),2)));
    y_0=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,-500<t_h&t_h<1),3),2)));
   
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,50<t_h&t_h<350),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_fef_vis_inout', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p1],3)'))

    p2=signrank(y_0,y_h)
    n=length(x_h);    
    eval(strcat('m_fef_vis_baseline', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(y_0), nanmean(y_h)-nanmean(y_0),nanstd(y_h-y_0)/sqrt(n),n, p2],3)'))

end

%%% IT Pref & Pref Correct
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_it_cr(ind_ses_it,1,:,:));
    var_cr4=squeeze(Power_it_cr(ind_ses_it,3,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,50<t_h&t_h<350),3),2)));
    y_0=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,-500<t_h&t_h<1),3),2)));
   
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,50<t_h&t_h<350),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_it_vis_pref', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n) ,n, p1],3)'))

    p2=signrank(y_0,y_h)
    n=length(x_h);    
    eval(strcat('m_it_vis_baseline', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(y_0), nanmean(y_h)-nanmean(y_0),nanstd(y_h-y_0)/sqrt(n) ,n, p2],3)'))

end

%% Late delay response stat test
%%% FEF IN & OUT Correct
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_fef_cr(ind_ses_fef,1,:,:));
    var_cr4=squeeze(Power_fef_cr(ind_ses_fef,2,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    y_0=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,-500<t_h&t_h<1),3),2)));
   
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_fef_del_inout', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p1],3)'))

    p2=signrank(y_0,y_h)
    n=length(x_h);    
    eval(strcat('m_fef_del_baseline', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(y_0), nanmean(y_h)-nanmean(y_0),nanstd(y_h-y_0)/sqrt(n),n, p2],3)'))

end

%%% IT Pref & Pref Correct
for ff=1:size(freq_bands,2)
    
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_cr1=squeeze(Power_it_cr(ind_ses_it,1,:,:));
    var_cr4=squeeze(Power_it_cr(ind_ses_it,3,:,:));
    
    x_h=func(squeeze(nanmean(nanmean(var_cr4(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    y_0=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,-500<t_h&t_h<1),3),2)));
   
    y_h=func(squeeze(nanmean(nanmean(var_cr1(:,f>f1&f<f2,t1<t_h&t_h<t2),3),2)));
    
    p1=signrank(x_h,y_h)
    n=length(x_h);    
    eval(strcat('m_it_del_pref', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n) ,n, p1],3)'))

    p2=signrank(y_0,y_h)
    n=length(x_h);    
    eval(strcat('m_it_del_baseline', bands_h_name{ff},'= round([nanmean(y_h) , nanmean(y_0), nanmean(y_h)-nanmean(y_0),nanstd(y_h-y_0)/sqrt(n) ,n, p2],3)'))

end

