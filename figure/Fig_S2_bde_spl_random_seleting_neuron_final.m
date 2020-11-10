%% Set Paths
clear all
clc
close all

%% Settings
load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')
load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
M1=0;
M2=0;

normalize_= @(x) x/nansum(x);
plv=@(x) abs(nanmean(exp(1j*x)));
% freqs2use  =[4:1:70 71:5:130];
% freqs2use  = 4:2:60;
%freqs2use= [3 4 10 12 15];

Color_loc=[1,0,0;0,0,1];
Color_obj=[1,0.5,0;0.5,0.5,0.5;0,1,0.5];
Color_bhv=[1,1,0;0,1,1];
win_=1;
t1=1100;t2=1700;N_bts=1001;

noramlized=0;


% freq_bands=[3.9,7.9,19.9,34.9,49.9;8.1,15.1,28.1,50.1,130.1];

freq_bands=[3.9,7.9,21.9,34.9,49.9,49.9;8.1,15.1,28.1,50.1,130.1,150.1];

%% load
% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\si_pval_it.mat')
load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\spl_ppc_per_session_bahmani_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\spl_ppc_per_session_bahmani_86_session.mat')

%% Cal INFO for neurons
spl_it_it_cr=[];spl_it_it_wr=[];
spl_fef_it_cr=[];spl_fef_it_wr=[];
spl_diff_it_cr=[];spl_diff_it_wr=[];
spl_it_it_n_cr=[];spl_it_it_n_wr=[];
spl_fef_it_n_cr=[];spl_fef_it_n_wr=[];
spl_diff_it_n_cr=[];spl_diff_it_n_wr=[];
pval_dsi_it=[];
pp_it=0;
for ss=1:91
    
    lfp=gen_lfp_pairs(ss);
    %     trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
    
    
    
    IT=[]; IT=resp(ss).IT;
    for nn=1:size(IT,2)
        obj_ix_h=[resp(ss).preff_obj(nn) resp(ss).npreff_obj(nn)];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        
        pp_it=pp_it+1;
        
        spl_it_it_cr(pp_it,:,:)=Info(ss).val_it_it_cr(nn,:,:);
        spl_it_it_wr(pp_it,:,:)=Info(ss).val_it_it_wr(nn,:,:);
        spl_fef_it_cr(pp_it,:,:)=Info(ss).val_fef_it_cr(nn,:,:);
        spl_fef_it_wr(pp_it,:,:)=Info(ss).val_fef_it_wr(nn,:,:);
        %       spl_diff_it_cr(pp_it,:,:)=Info(ss).val_diff_it_cr(nn,:,:);
        %       spl_diff_it_wr(pp_it,:,:)=Info(ss).val_diff_it_wr(nn,:,:);
        %
        spl_it_it_n_cr(pp_it,:,:)=Info(ss).val_it_it_cr(nn,:,:)-Info(ss).val_shuff_it_cr(nn,:,:);
        spl_it_it_n_wr(pp_it,:,:)=Info(ss).val_it_it_wr(nn,:,:)-Info(ss).val_shuff_it_wr(nn,:,:);
        %
        
        spl_fef_it_n_cr(pp_it,:,:)=Info(ss).val_fef_it_cr(nn,:,:)-Info(ss).val_shuff_it_cr(nn,:,:);
        spl_fef_it_n_wr(pp_it,:,:)=Info(ss).val_fef_it_wr(nn,:,:)-Info(ss).val_shuff_it_wr(nn,:,:);
        %       spl_diff_it_n_cr(pp_it,:,:)=Info(ss).val_diff_it_n_cr(nn,:,:);
        %       spl_diff_it_n_wr(pp_it,:,:)=Info(ss).val_diff_it_n_wr(nn,:,:);
        %
        ind_ses_it2(pp_it)=ss;
        
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition{:},obj_ix_h(1)+[0 3]),:);
        p_val_it(pp_it,1)=signrank(mean(var_h(:,t1:t2),2),mean(var_h(:,1:500),2));
        
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition{:},obj_ix_h(2)+[0 3]),:);
        p_val_it(pp_it,2)=signrank(mean(var_h(:,t1:t2),2),mean(var_h(:,1:500),2));
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition{:},obj_ix_h(3)+[0 3]),:);
        p_val_it(pp_it,3)=signrank(mean(var_h(:,t1:t2),2),mean(var_h(:,1:500),2));
        
        %            sig=[];  sig_=[];
        %             sig=resp(ss).IT{nn}(:,1:2100);
        %             N_trials=size(sig,1);
        %             for tt=1:size(sig,1)
        %                 sig_(tt,:)=SmoothN(sig(tt,:)*1, 10);
        %             end
        %
        %
        %          N_trials=size(resp(ss).IT{1},1);
        %           fr_in= mean2(sig_(logical(ismember(resp(ss).condition{:},[obj_ix_h(1) obj_ix_h(1)+3])),t1:t2));
        %             fr_out= mean2(sig_(logical(ismember(resp(ss).condition{:},[obj_ix_h(3) obj_ix_h(3)+3])),t1:t2));
        %             si=(fr_in-fr_out)/(fr_in+fr_out);
        %             for iteri=1:N_bts
        %                 ind_h=randperm(N_trials,round(N_trials/3));
        %                 fr_in_h= mean2(sig_(ind_h,t1:t2));
        %
        %                 ind_h_2=setdiff([1:N_trials],ind_h);
        %                 ind_h_ind=randperm(length(ind_h_2),round(N_trials/3));
        %                 ind_h=ind_h_2(ind_h_ind);
        %                 fr_out_h= mean2(sig_(ind_h,t1:t2));
        %                 si_h(iteri)=(fr_in_h-fr_out_h)/(fr_in_h+fr_out_h);
        %             end
        %             pval_dsi_it(pp_it)=1-sum(si>si_h)/N_bts;
        %
        
        
        
    end
    
end

spl_it_fef_cr=[];spl_it_fef_wr=[];
spl_fef_fef_cr=[];spl_fef_fef_wr=[];
spl_diff_fef_cr=[];spl_diff_fef_wr=[];

spl_it_fef_n_cr=[];spl_it_fef_n_wr=[];
spl_fef_fef_n_cr=[];spl_fef_fef_n_wr=[];
spl_diff_fef_n_cr=[];spl_diff_fef_n_wr=[];

pp_fef=0;
for ss=1:92
    
    lfp=gen_lfp_pairs(ss);
    %     trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
    
    
    obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
    obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
    
    FEF=[]; FEF=resp(ss).FEF;
    for nn=1:size(FEF,2)
        
        pp_fef=pp_fef+1;
        
        spl_it_fef_cr(pp_fef,:,:)=Info(ss).val_it_fef_cr(nn,:,:);
        spl_it_fef_wr(pp_fef,:,:)=Info(ss).val_it_fef_wr(nn,:,:);
        spl_fef_fef_cr(pp_fef,:,:)=Info(ss).val_fef_fef_cr(nn,:,:);
        spl_fef_fef_wr(pp_fef,:,:)=Info(ss).val_fef_fef_wr(nn,:,:);
        %       spl_diff_fef_cr(pp_fef,:,:)=Info(ss).val_diff_fef_cr(nn,:,:);
        %       spl_diff_fef_wr(pp_fef,:,:)=Info(ss).val_diff_fef_wr(nn,:,:);
        
        spl_it_fef_n_cr(pp_fef,:,:)=Info(ss).val_it_fef_cr(nn,:,:)-Info(ss).val_shuff_fef_cr(nn,:,:);
        spl_it_fef_n_wr(pp_fef,:,:)=Info(ss).val_it_fef_wr(nn,:,:)-Info(ss).val_shuff_fef_wr(nn,:,:);
        spl_fef_fef_n_cr(pp_fef,:,:)=Info(ss).val_fef_fef_cr(nn,:,:)-Info(ss).val_shuff_fef_cr(nn,:,:);
        spl_fef_fef_n_wr(pp_fef,:,:)=Info(ss).val_fef_fef_wr(nn,:,:)-Info(ss).val_shuff_fef_wr(nn,:,:);
        %       spl_diff_fef_n_cr(pp_fef,:,:)=Info(ss).val_diff_fef_n_cr(nn,:,:);
        %       spl_diff_fef_n_wr(pp_fef,:,:)=Info(ss).val_diff_fef_n_wr(nn,:,:);
        ind_ses_fef2(pp_fef)=ss;
        
    end
    
end

%% Session selection for comparing correct and Wrong
load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_full.mat');

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

control=1;

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


% %%%% Control for selecting one unit from each session
% if control
% idx_control=zeros(length(ind_ses_it2),1);
% ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_itphase)=1;
%
% for ss=1:92
%   ind_h= find(ismember(ind_ses_it2,ss)'&ind_resp);
%   if ~isempty(ind_h)
%   idx_control(ind_h(randi(length(ind_h),1)))=1;
%   end
% end
% ixx_it_itphase=find(ind_resp&idx_control);
% end
it_itphase_ind_m1=(ind_ses_it2(ixx_it_itphase)<52);
it_itphase_ind_m2=(ind_ses_it2(ixx_it_itphase)>51);



%

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


% ixx_it_=find(p_val_it(:,1)<0.05);
%  ixx_it=find(pval_dsi_it(idx_it_rate)<0.05);

% th_h=0.002;p_val_d=0.05;
%   ixx_it=idx_it;
ixx_fef=idx_fef;

ih_s=ixx_it(find((mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find((mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it_fefphase=setdiff(ixx_it,ih_s);

%%%% Control for selecting one unit from each session
ind_resp=[];ind_resp(ixx_it_fefphase)=1;

ses_index=ind_ses_it2(logical(ind_resp));
sesions_num= length(unique(ses_index));
unit_num=length(ses_index);


spl_fef_it_cr_tot=spl_fef_it_cr;
spl_fef_it_wr_tot=spl_fef_it_wr;
spl_fef_it_cr=[];spl_fef_it_cr=spl_fef_it_cr_tot(logical(ind_resp),:,:);
spl_fef_it_wr=[];spl_fef_it_wr=spl_fef_it_wr_tot(logical(ind_resp),:,:);
%%
for iteri=1:1000
    
    if control
        idx_control=zeros(length(ses_index),1);
        %  ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_fefphase)=1;
        num_=1;
        for ss=1:92
            ind_h=[];ind_h= find(ismember(ses_index,ss));
            if ~isempty(ind_h)
                if length(ind_h)>=num_
                    idx_control(ind_h(randperm(length(ind_h),randi(num_))))=1;
                else
                    idx_control(ind_h(randperm(length(ind_h),length(ind_h))))=1;
                end
                
            end
        end
        %         if sum(idx_control)<134
        %         non_sel=find(~idx_control);
        %         n_h=134-sum(idx_control);
        %        idx_control( non_sel(randperm(length(non_sel),n_h)))=1;
        %         end
        ixx_it_fefphase=find(idx_control);
        %         (sum(idx_control)-(((sum(ind_resp)/82)-1)*82))
    end
    it_fefphase_ind_m1=(ses_index(ixx_it_fefphase)<52);
    it_fefphase_ind_m2=(ses_index(ixx_it_fefphase)>51);
    
    
    %
    noramlized=0;
    
    
    %% Scatter Cr versus Wr In preff
    in_diff=[];
    condi=1;f=freqs2use;
    % delta_=(bb2_-bb1_)/10;
    %% Scatter Correct versus  Wrong
    bin_=0.5;
    hist_var=[];nbin1=-2:bin_:2.5;nbin2=-2.5:bin_:2;
    ff=3;
    f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
    
    var_wr = [];var_wr = squeeze((spl_fef_it_wr(ixx_it_fefphase,1,:)));
    var_cr = [];var_cr = squeeze((spl_fef_it_cr(ixx_it_fefphase,1,:)));
    
    x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));
    
    y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2)));
    
    in_diff(ff,:)=y_h-x_h;
    p=signrank(x_h,y_h);
    mean_val(iteri)=nanmedian(y_h)-nanmedian(x_h);
    p_val(iteri)=p;
    numb(iteri)=length(x_h);
end

%% Anova
val_h_d=spl_fef_it_cr-spl_fef_it_wr;
v_h=(squeeze((nanmean(val_h_d(:,1,f>f1&f<f2),3))));
same_pairs=[];pp=0;
for ss=1:92
    ind_h=[];ind_h= find(ismember(ses_index,ss));
    if length(ind_h)>1
        P_h= nchoosek(ind_h,2);
        if size(P_h,1)>0
            if size(P_h,1)<2
                
                ll= size(same_pairs,1);
                same_pairs(ll+1:ll+1,:)= [[v_h(P_h)]];
            else
                
                same_pairs= [same_pairs; v_h(P_h)];
                
            end
        end
    end
end
%%
figure('name','Control for Fig. 4a')
hold on
scatter(same_pairs(:,1),same_pairs(:,2))
export_file_to_excel(same_pairs,'figS2d',10)

x_h=same_pairs(:,1);y_h=same_pairs(:,2);
x_h_=x_h(~(isnan(y_h)&isnan(x_h)));
y_h_=y_h(~(isnan(y_h)&isnan(x_h)));
x_h=x_h_;y_h=y_h_;
% [r,p]=correlation(x_h,y_h,1000)
x=x_h;y=y_h;
tbl = table(x,y);
lm = fitlm(tbl,'linear');
b=lm.Coefficients{1,1};
m=lm.Coefficients{2,1};
xx=-0.15:0.1:.15;
% [r,p]=corr(x_h,y_h,'type','kendal');
% [r,p]=corr(x_h,y_h,'type','kendal')
[r,p]=correlation(x_h,y_h,10)
cohen_r2=(r^2)

length(x_h)

plot(xx,m*xx+b,'r')


xlabel(' SPL Cr-Wr[neuron1](a.u.)')
ylabel(' SPL Cr-Wr[neuron2](a.u.)')
ylim([-0.16 0.16])
xlim([-0.16 0.16])
% line([0.5 0.5],ylim)
% line(xlim,[0 0])
set(gca,'fontsize',14,'fontweight','bold')
text( 0.10,0.12,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.10,0.10,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.10,0.08,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')

%%

P_h= nchoosek(ses_index,2);
r_a=[];
for iteri=1:100
    ind_h=randperm(size(P_h,1),234);
    val_h= v_h(P_h(ind_h,:));
    x_h=val_h(:,1); y_h=val_h(:,2);
    x_h_=x_h(~(isnan(y_h)&isnan(x_h)));
    y_h_=y_h(~(isnan(y_h)&isnan(x_h)));
    x_h=x_h_;y_h=y_h_;
    [r_a(iteri),p]=correlation(x_h,y_h,10);
    
end
%% random pairs 
figure
hold on
[rr cc]=hist(r_a,-1:0.05:1);
export_file_to_excel([cc',[rr/sum(rr)]'],'figS2e',9)

bar(cc,rr/sum(rr),'barwidth',1)
% [p,h,stats] = signrank(r_a,r)
line([r r],ylim,'color','r','linewidth',2)
set(gca,'fontsize',14,'fontweight','bold')
xlim([-0.25 0.25])
xlabel('Correlation value (a.u.)')
ylabel('% pairs')
text( r,0.3,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( r,0.25,strcat('p = ',num2str(sum(r_a>r)/length(r_a),2)),'fontsize',14,'fontweight','bold')

%%
figure('name','Control for Fig. 4a')
bins=0.0:0.005:0.5;
[counts,centers]=hist(p_val,bins);
counts=counts/sum(counts)*100;
bar(centers,counts,1)
line([0 0],ylim)
line([0.05 0.05],ylim)
set(gca,'fontsize',14,'fontweight','bold')
text( 0.05,5,strcat(num2str((sum(p_val<0.05)/length(p_val)*100)),'% of valuse <0.05'))

ylabel('%')
xlabel('p values (a.u.)')
xlim([0 0.5])
title('Control for Fig. 4a (200 times)')


figure('name','fig S2b Control for Fig. 4a')
bins=0.0:0.005:0.5;
[counts,centers]=hist(mean_val);
counts=counts/sum(counts)*100;
export_file_to_excel([centers',counts'],'figS2b',9)

bar(centers,counts,1)
line([0 0],ylim)
% line([0.05 0.05],ylim)
signrank(mean_val)
line([median(mean_val) median(mean_val)],ylim)
set(gca,'fontsize',14,'fontweight','bold')
text( 0.0,5,strcat('p=',num2str(signrank(mean_val,0)),'% of valuse <0.05'))
text( 0.05,5,strcat(num2str((sum(mean_val>0.0)/length(mean_val)*100))))

ylabel('%')
xlabel('Cr - Wr (a.u.)')
xlim([-0.01 0.03])
title('Control for Fig. 4a ')


%%  correlation between
cal=0
if cal
    for pp=1:1000
        ind_s_pp=randperm(92,60);
        mean_val_rand=[];mean_val=[];
        for iteri=1:100
            
            idx_control=zeros(length(ses_index),1);
            %  ses=[];ind_resp=zeros(253,1);ind_resp(ixx_it_fefphase)=1;
            num_=1;
            for ss=ind_s_pp
                ind_h=[];ind_h= find(ismember(ses_index,ss));
                if ~isempty(ind_h)
                    if length(ind_h)>=num_
                        idx_control(ind_h(randperm(length(ind_h),randi(num_))))=1;
                    else
                        idx_control(ind_h(randperm(length(ind_h),length(ind_h))))=1;
                    end
                    
                end
            end
            %         if sum(idx_control)<134
            %         non_sel=find(~idx_control);
            %         n_h=134-sum(idx_control);
            %        idx_control( non_sel(randperm(length(non_sel),n_h)))=1;
            %         end
            ixx_it_fefphase=find(idx_control);
            
            ixx_it_fefphase_rand=find(ismember(ses_index,ind_s_pp));
            
            
            %
            noramlized=0;
            
            
            %%Scatter Cr versus Wr In preff
            in_diff=[];
            condi=1;f=freqs2use;
            % delta_=(bb2_-bb1_)/10;
            %%Scatter Correct versus  Wrong
            bin_=0.5;
            hist_var=[];nbin1=-2:bin_:2.5;nbin2=-2.5:bin_:2;
            ff=3;
            f1=freq_bands(1,ff);  f2=freq_bands(2,ff);
            
            var_wr = [];var_wr = squeeze((spl_fef_it_wr(ixx_it_fefphase,1,:)));
            var_cr = [];var_cr = squeeze((spl_fef_it_cr(ixx_it_fefphase,1,:)));
            x_h=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));
            y_h=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2)));
            %     in_diff(ff,:)=y_h-x_h;
            p=signrank(x_h,y_h);
            mean_val(iteri)=nanmedian(y_h)-nanmedian(x_h);
            mean_val_d(iteri)=(nanmean(y_h)-nanmean(x_h))*sqrt(length(y_h))/nanstd(y_h-x_h);
            
            p_val(iteri)=p;
            numb(iteri)=length(x_h);
            
            
            var_wr = [];var_wr = squeeze((spl_fef_it_wr(ixx_it_fefphase_rand,1,:)));
            var_cr = [];var_cr = squeeze((spl_fef_it_cr(ixx_it_fefphase_rand,1,:)));
            x_h_rand=(squeeze(nanmean(nanmean(var_wr(:,f>f1&f<f2),2),2)));
            y_h_rand=(squeeze(nanmean(nanmean(var_cr(:,f>f1&f<f2),2),2)));
            %     in_diff_rand(ff,:)=y_h_rand-x_h_rand;
            p=signrank(x_h_rand,y_h_rand);
            mean_val_rand(iteri)=nanmedian(y_h_rand)-nanmedian(x_h_rand);
            mean_val_d_rand(iteri)=(nanmean(y_h_rand)-nanmean(x_h_rand))*sqrt(length(y_h_rand))/nanstd(y_h_rand-x_h_rand);
            
            p_val_rand(iteri)=p;
            
        end
        val1(pp)=nanmean(mean_val_rand);
        val2(pp)=nanmean(mean_val);
        val1_d(pp)=nanmean(mean_val_d_rand);
        val2_d(pp)=nanmean(mean_val_d);
        
    end
    %%
    figure
    scatter(val1_d,val2_d)
    [r,p]=corr(val1_d',val2_d')
    text( 1.5,1,strcat('p=',num2str(p)),'fontsize',14,'fontweight','bold')
    text( 1.5,1.5,strcat('r=',num2str(r)),'fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
    % axis([0 0.02 0 0.02])
    ylabel('Independent')
    xlabel('All ')
    % xlim([-0.03 0.03])
    title('Control for Fig. 4a ')
    
    
    
    figure
    scatter(val1,val2)
    [r,p]=corr(val1',val2')
    text( 0.015,0.01,strcat('p=',num2str(p)),'fontsize',14,'fontweight','bold')
    text( 0.015,0.015,strcat('r=',num2str(r)),'fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
    % axis([0 0.02 0 0.02])
    ylabel('Independent')
    xlabel('All ')
    % xlim([-0.03 0.03])
    title('Control for Fig. 4a ')
    
end
