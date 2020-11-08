function Fig5_Object_coding(resp)

%% In the name of Allah
% Object coding in IT

%% Set Paths
% clear all
% close all
% clc

%% Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')

%% Settings
ind_b=[1 500];
Ns=size(resp,2);
win_=300;step=10;
Percenti=33.33;
control=0;
Shuffle_highlow_PPL=0;
Shuffle_Preferance=1;
baseline_corrected_plv=1;
N_shuff=1001;
%ROC_,  F_ratio_base, M_index mi_object_coding
%      inf_cal=@(x,y,z) F_ratio_base(x,y,z);  method='vep';tresh=-10;  inf_pre=@(x,y) (x-y)./(y+x+1);
%      inf_cal=@(x,y,z) mi_object_coding(x,y,z); method='mi';
       inf_cal=@(x,y,z) ROC_(x,y); method='roc';inf_pre=@(x,y) (x-y)./(y+x);
%      inf_cal=@(x,y,z) ROC_bts(x,y); method='roc_bts';
%      inf_cal=@(x,y,z) (nanmean(x)-nanmean(y))/sqrt(((nanstd(x)^2)/length(x))+((nanstd(y)^2)/length(y)));method='dprim';
%      inf_cal=@(x,y,z) abs(nanmean(x)-nanmean(y))/sqrt(((nanstd(x)^2)/length(x))+((nanstd(y)^2)/length(y))); method='abs_dprim';
  %ROC_  F_ratio_base
    %ROC_  F_ratio_base
  
%% Calculate PPL single trials
cal=1;
if cal
    f1=19;f2=32; bw=3;%   f1=18;f2=28;
    t1_ld=600;t2_ld=1200;
    st_ind=1;end_ind=2100;
    t1_b=-500;t2_b=1;
    load (strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\plv_wavelet_over_time_full.mat'))
    freqs2use_plv = freqs2use;
    shuffle_corrected = 0;
    if shuffle_corrected
        plv_cr = (plv_cr-plv_shuff_m)./plv_shuff_std;plv_wr = (plv_wr-plv_shuff_m)./plv_shuff_std;
        %     plv_cr = (plv_shuff_m);plv_wr = plv_shuff_std;
    end
    %%Baseline Correction
    baseline_corrected = 0;
    if baseline_corrected
        ind_b2 = [1 30];
        plv_cr_ = plv_cr;plv_wr_ = plv_wr;
        val_b = nanmean(plv_cr_(:,:,:,ind_b2(1):ind_b2(2)),4);
        val_sd = nanstd(val_b);
        val_sd = repmat(val_sd,NS,1,1);
        val_sd = repmat(val_sd,1,1,1,size(plv_cr,4));
        val_b = repmat(val_b,1,1,1,size(plv_cr,4));
        plv_cr = (plv_cr_-val_b)./val_sd;
        
        val_b = nanmean(plv_wr_(:,:,:,ind_b(1):ind_b(2)),4);
        val_sd = nanstd(val_b);
        val_sd = repmat(val_sd,NS,1,1);
        val_sd = repmat(val_sd,1,1,1,size(plv_cr,4));
        val_b = repmat(val_b,1,1,1,size(plv_cr,4));
        plv_wr = (plv_wr_-val_b)./val_sd;
    end
    f=freqs2use_plv;
    plv_all=[];plv_all=squeeze((nanmean(plv_cr(:,:,:,t1_ld<t_h&t_h<t2_ld),4)));
    plv_cat_wav=[];plv_cat_hilb=[];
    pp_it=0;pp_fef=0;it_tag=[];
    for ss=1:Ns
        tic
        
        %   plv_val=plv_all(:,1:3,freqs2use_plv>f1&freqs2use_plv<f2);
        plv_val=plv_all(:,:,freqs2use_plv>f1&freqs2use_plv<f2);
        
        [val_max,ind_f]= max(nanstd(squeeze(plv_val(ss,:,:)),[],1));
        fw1=f1+ind_f-bw;  fw2=f1+ind_f+bw;
        if ind_f<2 fw1=f1+1;fw2=f1+bw-2;end
        if ind_f>size(plv_val,3)-1; fw1=size(plv_val,3)+f1-bw-2;fw2=size(plv_val,3)+f1; end
        
        %         [val_max,ind_f]= max(plv_val(ss,:));ind_max=ind_f+fw1-1;
        %        [~ ,in_bw1]= min(plv_val(ss,1:ind_f)-(val_max*0.85));
        %        [~ ,in_bw2]= min(plv_val(ss,ind_f:end)-(val_max*0.85));
        %
        %
        %       f1=in_bw1+fw1-1; f2=ind_f+in_bw2+fw1-1;
        % if ind_max==fw2;f1=fw2-4;f2=fw2+4;end
        % if ind_max==fw1;f1=fw1-4;f2=fw2-4;end
        
        %             f1=ind_f-BW;f2=ind_f+BW; %   f1=18;f2=28;
        %         if ss<52
        %         fw1=22;fw2=28;
        %         else
        %         fw1=20;fw2=27;
        %
        %         end
        fw1=22;fw2=28;
        freqs2use= fw1:fw2;
        ss
        lfp=gen_lfp_pairs(ss);
        trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
        s_it = []; s_it = lfp.IT(trii_h,st_ind:end_ind);
        s_fef = []; s_fef = lfp.FEF(trii_h,st_ind:end_ind);
        %% Wavelet
        sig_it = [];sig_fef = [];
        [sig_it,sig_fef,freqs2use,t]= wavelet_core_pertrial(s_it,s_fef,freqs2use);
        t_h=t-500;trial_num=size(sig_it,1);
        phi_fef=[];phi_fef=angle(sig_fef);
        phi_it=[];phi_it=angle(sig_it);
        diff_phi=angle(exp(1j*(phi_fef-phi_it)));
        
        plv_baseline=[];
        plv_baseline=squeeze(abs(mean((exp(1j*(diff_phi(:,:,t1_b<t_h&t_h<t2_b)))),3)));
        plv_delay=[];
        plv_delay=squeeze(abs(mean((exp(1j*(diff_phi(:,:,t1_ld<t_h&t_h<t2_ld)))),3)));
        %%% Normalized with baseline
        %         %         plv_h_map= squeeze(val_(:,2,:)-repmat(nanmean(val_(:,1,:)),size(val_,1),1,1));
        %         plv_h_map= squeeze(val_(:,2,:));
        %         plv_h =[]; plv_h= nanmean(plv_h_map(:,freqs2use>=f1&freqs2use<=f2),2);
        
        win_h=[];
        time_window_idx = 100;%round((1000/22)*2);
        for ti=1+time_window_idx:10:1700-time_window_idx;win_h=[win_h; [ti-time_window_idx ti+time_window_idx]];end
        tt_h=mean(win_h ,2);tt_h=tt_h-500;
        plv_h_all=[];
        for ti = 1:size(win_h,1)
            var_h_w= [];var_h_w=squeeze(abs(nanmean(exp(1j*(phi_fef(:,:,win_h(ti,1):win_h(ti,2))-...
                phi_it(:,:,win_h(ti,1):win_h(ti,2)))),3)));
            plv_h_all(:,:,ti)=var_h_w;
            
        end
        
        %         val_b=(nanmean(plv_h_cr(:,:,1:40),3));
        %
        %         val_b_m=nanmean(val_b);val_b_m=repmat(val_b_m,150,1)';
        %         val_b_m_(1,:,:)=val_b_m;val_b_m=repmat(val_b_m_,trial_num,1,1);
        %         val_b_std=nanmean(val_b);val_b_std=repmat(val_b_std,150,1)';
        %         val_b_std_(1,:,:)=val_b_std;val_b_std=repmat(val_b_std_,trial_num,1,1);
        %         plv_h_cr=(plv_h_cr-val_b_m)./val_b_std;
        
        %         plv_h =[]; plv_h= squeeze((nanmean(plv_h_cr(:,:,(600+time_window_idx)<tt_h&tt_h<(1200-time_window_idx)),3)));
        %         plv_std=[];
        %               plv_std=[];for ff=2:size(plv_h,2)-1; plv_std(ff)=nanstd(reshape(plv_h(:,ff-1:ff+1),[],1));end
        
        %         for ff=2:size(plv_h,2)-1;
        %             [p, tb1]=anovan(nanmean(plv_h(:,ff-1:ff+1),2),{lfp.condition{:}(trii_h)},'display','off');
        %             plv_std(ff)=tb1{2,7};
        %         end
        %         plv_std=plv_std(2:end);
        %
        %
        %         [val_max,ind_f]= max(plv_std);ind_max=ind_f+fw1-1;
        %         in_bw1=ind_f-4;        in_bw2=ind_f+4;
        %
        %         [~ ,in_bw1]= min(plv_std(1:ind_f)-(val_max*0.85));
        %         [~ ,in_bw2]= min(plv_std(ind_f:end)-(val_max*0.85));
        %
        %
        %         if ind_f>=length(plv_std)-3;in_bw1=length(plv_std)-4;in_bw2=length(plv_std);end
        %         if ind_f<3;in_bw1=1;in_bw2=5;end
        %         if in_bw1<1;in_bw1=1;end
        %         if in_bw2>length(plv_std);in_bw2=length(plv_std);end
        %
        
        %         plv_h =[]; plv_h= squeeze(nanmean(nanmean(plv_h_cr(:,in_bw1:in_bw2,(600+time_window_idx)<tt_h&tt_h<(1200-time_window_idx)),3),2));
        plv_h =[]; plv_h= squeeze(nanmean(nanmean(plv_h_all(:,:,(600+time_window_idx)<tt_h&tt_h<(1200-time_window_idx)),3),2));
        t=1:2100;t=t-500;
        pow_it_delay =[]; pow_it_delay= squeeze((nanmean(abs(sig_it(:,:,(t1_ld)<t&t<(t2_ld)).^2),3)));
        pow_fef_delay =[]; pow_fef_delay= squeeze((nanmean(abs(sig_fef(:,:,(t1_ld)<t&t<(t2_ld)).^2),3)));
        
        pow_it_baseline =[]; pow_it_baseline= squeeze((nanmean(abs(sig_it(:,:,(t1_b)<t&t<(t2_b)).^2),3)));
        pow_fef_baseline =[]; pow_fef_baseline= squeeze((nanmean(abs(sig_fef(:,:,(t1_b)<t&t<(t2_b)).^2),3)));
       
        plv_b =[]; plv_b= squeeze((nanmean(plv_h_all(:,:,tt_h<0),3)));
        plv_hd =[]; plv_hd= squeeze((nanmean(plv_h_all(:,:,(600+time_window_idx)<tt_h&tt_h<(1200-time_window_idx)),3)));
        %         plv_cat_wav{ss}=nanmean(plv_hd-plv_b,2);
        %         plv_cat_wav{ss}=nanmean(plv_delay-plv_baseline,2);
        
        if baseline_corrected_plv
        plv_cat_wav{ss}=nanmean(plv_delay-plv_baseline,2);
        pow_it_wav{ss}=nanmean(pow_it_delay-pow_it_baseline,2);
        pow_fef_wav{ss}=nanmean(pow_fef_delay-pow_fef_baseline,2);
        else
        plv_cat_wav{ss}=nanmean(plv_delay,2);
        pow_it_wav{ss}=nanmean(pow_it_delay,2);
        pow_fef_wav{ss}=nanmean(pow_fef_delay,2);
        end
        freqs_of_ses{ss}=freqs2use;
    end
%     save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_wideband_full_test.mat',...
%         'plv_cat_hilb','plv_cat_wav','pow_it_wav','pow_fef_wav','freqs_of_ses')
    
end

%% Rate categorized calculate
cal=1;
if cal
%     load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_wideband_full_test.mat')
    %         load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_gammaband_full.mat')
    
    % load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_20_28_full.mat')
    t=1:2100;
    t_end=2100;
    
    time_win=[];
    %     time_win=[1100 1700];
    for tt=1+win_/2:step:2100-win_/2;time_win =[time_win; [tt-win_/2 tt+win_/2]];end
    %     time_win=time_win(time_win(:,2)<2100,:);
    pref_chang_both=zeros(200,1);
    pref_chang_in=zeros(200,1);
    pref_chang_out=zeros(200,1);
    
    for mpi=1:3
        mi_all_in=[]; mi_high_in=[];mi_low_in=[];
        mi_all_out=[]; mi_high_out=[];mi_low_out=[];
        mi_all_both=[]; mi_high_both=[];mi_low_both=[];
        
        rate_all_in=[]; rate_high_in=[];rate_low_in=[];
        rate_all_out=[]; rate_high_out=[];rate_low_out=[];
        rate_all_both=[]; rate_high_both=[];rate_low_both=[];
        selec_all=1;
        
        pp_it=0;pp_fef=0;it_tag=[];
        for ss=1:Ns
            lfp=gen_lfp_pairs(ss);
            trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            
            
            
            t_h_sp=t-500;
            IT=resp(ss).IT;
            if ~isempty(IT)
                for nn=1:size(IT,2)
                    obj_ix_h=[resp(ss).preff_obj(nn) resp(ss).npreff_obj(nn)];
                    obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
                    it_tag=[it_tag ss*10+nn];
                    pp_it=pp_it+1;
                    it_=IT{nn}(trii_h,:);
                    
                    sig_=[];
                    for tt=1:size(it_,1)
                        sig_(tt,:)=SmoothN(it_(tt,:)*1, 10);
                    end
                    
                    var_b=(mean2(sig_(:,ind_b(1):ind_b(2))));
                    var_max=max(max(sig_));
                    sig_=(sig_-var_b)./(var_max-var_b);
                    
                    %  rate_= nanmean(sig_(:,t1_ld<t_h_sp&t_h_sp<t2_ld),2);
                    %      rate_= (sig_(:,t1_ld<t_h_sp&t_h_sp<t2_ld));
                    rate_= sig_;
                    
                    for condi=1:6
                        if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
                        
                        %Correct
                        ix_h = [];
                        ix_h = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==1)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                        
                        %                 plv_cat=plv_cat_hilb{ss};
                        plv_cat=[];
                        switch mpi
                            case 1
                                plv_cat=plv_cat_wav{ss};
                                
                            case 2
                                plv_cat=pow_it_wav{ss};
                                
                            case 3
                                plv_cat=pow_fef_wav{ss};
                        end
                        plv_=plv_cat(ix_h);
                        rate_h_=rate_(ix_h,:);
                        plv_high_h=(plv_>prctile(plv_,round(100-Percenti)));
                        plv_low_h=(plv_<prctile(plv_,round(Percenti)));
                        
                        
                        val_it_all{pp_it,condi}=rate_h_;
                        val_it_low{pp_it,condi}=rate_h_(plv_low_h,:);
                        val_it_high{pp_it,condi}=rate_h_(plv_high_h,:);
                        
                    end
                end
            end
            ss
        end
        
        
        t1_prefeance=1;t2_prefeance=1300;%    tt_1=1100;tt_2=1700;
        pref_chang=0;
        
        t_end=2100;%size(val_it_all{1,1},2);
        
        time_win=[];
        %     time_win=[1100 1700];
        for tt=1+win_/2:step:2100-win_/2;time_win =[time_win; [tt-win_/2 tt+win_/2]];end
        %     time_win=time_win(time_win(:,2)<2100,:);
        pref_chang_both=zeros(200,1);
        pref_chang_in=zeros(200,1);
        pref_chang_out=zeros(200,1);
        
        mi_all_in=[]; mi_high_in=[];mi_low_in=[];
        mi_all_out=[]; mi_high_out=[];mi_low_out=[];
        mi_all_both=[]; mi_high_both=[];mi_low_both=[];
        
        mi_all_in_n=[]; mi_high_in_n=[];mi_low_in_n=[];
        mi_all_out_n=[]; mi_high_out_n=[];mi_low_out_n=[];
        mi_all_both_n=[]; mi_high_both_n=[];mi_low_both_n=[];
        
        mi_all_in_p=[]; mi_high_in_p=[];mi_low_in_p=[];
        mi_all_out_p=[]; mi_high_out_p=[];mi_low_out_p=[];
        mi_all_both_p=[]; mi_high_both_p=[];mi_low_both_p=[];
        
        rate_all_in=[]; rate_high_in=[];rate_low_in=[];
        rate_all_out=[]; rate_high_out=[];rate_low_out=[];
        rate_all_both=[]; rate_high_both=[];rate_low_both=[];
        
        mi_all_in_shuff=[]; mi_high_in_shuff=[];mi_low_in_shuff=[];
        mi_all_out_shuff=[]; mi_high_out_shuff=[];mi_low_out_shuff=[];
        mi_all_both_shuff=[]; mi_high_both_shuff=[];mi_low_both_shuff=[];
        
        rate_all_in_shuff=[]; rate_high_in_shuff=[];rate_low_in_shuff=[];
        rate_all_out_shuff=[]; rate_high_out_shuff=[];rate_low_out_shuff=[];
        rate_all_both_shuff=[]; rate_high_both_shuff=[];rate_low_both_shuff=[];
        
        
        selec_all=1;
        for pp=1:size(val_it_all,1)
            pp
            tic
            
            tt=1;
            R1_p=[]; R2_p=[]; R3_p=[];
            %         R1_p= [nanmean(val_it_all{pp,1}(:,t1_prefeance:t2_prefeance),2);nanmean(val_it_all{pp,4}(:,t1_prefeance:t2_prefeance),2)];
            %         R2_p= [nanmean(val_it_all{pp,2}(:,t1_prefeance:t2_prefeance),2);nanmean(val_it_all{pp,5}(:,t1_prefeance:t2_prefeance),2)];
            %         R3_p= [nanmean(val_it_all{pp,3}(:,t1_prefeance:t2_prefeance),2);nanmean(val_it_all{pp,6}(:,t1_prefeance:t2_prefeance),2)];
            %
            R1_p= [nanmean(val_it_all{pp,1}(:,t1_prefeance:t2_prefeance),2)];
            R2_p= [nanmean(val_it_all{pp,2}(:,t1_prefeance:t2_prefeance),2)];
            R3_p= [nanmean(val_it_all{pp,3}(:,t1_prefeance:t2_prefeance),2)];
            
            if pref_chang
                ind_obj_max=[];ind_obj_min=[];
                
                %             [~,ind_obj_max]= max([nanmedian(R1_p) nanmedian(R2_p) nanmedian(R3_p)]);
                %             [~,ind_obj_min]= min([nanmedian(R1_p) nanmedian(R2_p) nanmedian(R3_p)]);
                %
                [~,ind_obj_max]= max([nanmean(R1_p) nanmean(R2_p) nanmean(R3_p)]);
                [~,ind_obj_min]= min([nanmean(R1_p) nanmean(R2_p) nanmean(R3_p)]);
                interm=setdiff([1:3],[ind_obj_max,ind_obj_min]);
                
            else
                ind_obj_max=1;ind_obj_min=3;interm=2;
            end
            
            
            N1_i=size(val_it_all{pp,ind_obj_max},1);
            N2_i=size(val_it_all{pp,interm},1);
            N3_i=size(val_it_all{pp,ind_obj_min},1);
            N1_o=size(val_it_all{pp,ind_obj_max+3},1);
            N2_o=size(val_it_all{pp,interm+3},1);
            N3_o=size(val_it_all{pp,ind_obj_min+3},1);
            ind1_i=[];ind2_i=[];ind3_i=[];
            ind1_o=[];ind2_o=[];ind3_o=[];
            ind1_b=[];ind2_b=[];ind3_b=[];
            for iteri=1:N_shuff
                ind1_i(iteri,:)=randperm(N1_i);
                ind2_i(iteri,:)=randperm(N2_i);
                ind3_i(iteri,:)=randperm(N3_i);
                ind1_o(iteri,:)=randperm(N1_o);
                ind2_o(iteri,:)=randperm(N2_o);
                ind3_o(iteri,:)=randperm(N3_o);
                ind1_b(iteri,:)=randperm(N1_i+N1_o);
                ind2_b(iteri,:)=randperm(N2_i+N2_o);
                ind3_b(iteri,:)=randperm(N3_i+N3_o);
            end
            
            for tt=1:size(time_win,1)
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_both(pp,tt)=mi;
                rate_all_both(pp,1,tt)=nanmean(R1);
                rate_all_both(pp,2,tt)=nanmean(R2);
                rate_all_both(pp,3,tt)=nanmean(R3);
                
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_both(pp,tt)=mi;
                rate_high_both(pp,1,tt)=nanmean(R1);
                rate_high_both(pp,2,tt)=nanmean(R2);
                rate_high_both(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_both(pp,tt)=mi;
                rate_low_both(pp,1,tt)=nanmean(R1);
                rate_low_both(pp,2,tt)=nanmean(R2);
                rate_low_both(pp,3,tt)=nanmean(R3);
                
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all = [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all = [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all = [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        ind1=ind1_b(iteri,:);
                        ind2=ind2_b(iteri,:);
                        ind3=ind3_b(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_both_shuff(pp,tt,iteri)=mi;
                        rate_low_both_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_both_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_both_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_both_shuff(pp,tt,iteri)=mi;
                        rate_high_both_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_both_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_both_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                %         end
                %
                %
                %         for tt=1:size(time_win,1)
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_in(pp,tt)=mi;
                rate_all_in(pp,1,tt)=nanmean(R1);
                rate_all_in(pp,2,tt)=nanmean(R2);
                rate_all_in(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_in(pp,tt)=mi;
                rate_high_in(pp,1,tt)=nanmean(R1);
                rate_high_in(pp,2,tt)=nanmean(R2);
                rate_high_in(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_in(pp,tt)=mi;
                rate_low_in(pp,1,tt)=nanmean(R1);
                rate_low_in(pp,2,tt)=nanmean(R2);
                rate_low_in(pp,3,tt)=nanmean(R3);
                
                
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        ind1=ind1_i(iteri,:);
                        ind2=ind2_i(iteri,:);
                        ind3=ind3_i(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_in_shuff(pp,tt,iteri)=mi;
                        rate_low_in_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_in_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_in_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_in_shuff(pp,tt,iteri)=mi;
                        rate_high_in_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_in_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_in_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                
                %         end
                %
                %
                %         for tt=1:size(time_win,1)
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_out(pp,tt)=mi;
                rate_all_out(pp,1,tt)=nanmean(R1);
                rate_all_out(pp,2,tt)=nanmean(R2);
                rate_all_out(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_out(pp,tt)=mi;
                rate_high_out(pp,1,tt)=nanmean(R1);
                rate_high_out(pp,2,tt)=nanmean(R2);
                rate_high_out(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_out(pp,tt)=mi;
                rate_low_out(pp,1,tt)=nanmean(R1);
                rate_low_out(pp,2,tt)=nanmean(R2);
                rate_low_out(pp,3,tt)=nanmean(R3);
                
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all= [nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all= [nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all= [nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        
                        ind1=ind1_o(iteri,:);
                        ind2=ind2_o(iteri,:);
                        ind3=ind3_o(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_out_shuff(pp,tt,iteri)=mi;
                        rate_low_out_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_out_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_out_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_out_shuff(pp,tt,iteri)=mi;
                        rate_high_out_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_out_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_out_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                
            end
            toc
        end
        t_h=mean(time_win,2)-500;
        
        switch mpi
            case 1
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat'),...
                    'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
                    'mi_low_in','mi_low_out','mi_low_both','t_h')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat'),...
                    'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
                    'rate_low_in','rate_low_out','rate_low_both','t_h')
            case 2
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat'),...
                    'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
                    'mi_low_in','mi_low_out','mi_low_both','t_h')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat'),...
                    'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
                    'rate_low_in','rate_low_out','rate_low_both','t_h')
            case 3
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat'),...
                    'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
                    'mi_low_in','mi_low_out','mi_low_both','t_h')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat'),...
                    'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
                    'rate_low_in','rate_low_out','rate_low_both','t_h')
        end
        if Shuffle_highlow_PPL
            
            
            save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_wideband_full_shuffled.mat',...
                'mi_high_in_shuff','mi_high_out_shuff','mi_high_both_shuff',...
                'mi_low_in_shuff','mi_low_out_shuff','mi_low_both_shuff',...
                't_h')
        end
        
        
        
        %%%% Interval
        
        time_win=[];
        time_win=[1100 1700];
        mi_all_in=[]; mi_high_in=[];mi_low_in=[];
        mi_all_out=[]; mi_high_out=[];mi_low_out=[];
        mi_all_both=[]; mi_high_both=[];mi_low_both=[];
        
        mi_all_in_n=[]; mi_high_in_n=[];mi_low_in_n=[];
        mi_all_out_n=[]; mi_high_out_n=[];mi_low_out_n=[];
        mi_all_both_n=[]; mi_high_both_n=[];mi_low_both_n=[];
        
        mi_all_in_p=[]; mi_high_in_p=[];mi_low_in_p=[];
        mi_all_out_p=[]; mi_high_out_p=[];mi_low_out_p=[];
        mi_all_both_p=[]; mi_high_both_p=[];mi_low_both_p=[];
        
        rate_all_in=[]; rate_high_in=[];rate_low_in=[];
        rate_all_out=[]; rate_high_out=[];rate_low_out=[];
        rate_all_both=[]; rate_high_both=[];rate_low_both=[];
        
        mi_all_in_shuff=[]; mi_high_in_shuff=[];mi_low_in_shuff=[];
        mi_all_out_shuff=[]; mi_high_out_shuff=[];mi_low_out_shuff=[];
        mi_all_both_shuff=[]; mi_high_both_shuff=[];mi_low_both_shuff=[];
        
        rate_all_in_shuff=[]; rate_high_in_shuff=[];rate_low_in_shuff=[];
        rate_all_out_shuff=[]; rate_high_out_shuff=[];rate_low_out_shuff=[];
        rate_all_both_shuff=[]; rate_high_both_shuff=[];rate_low_both_shuff=[];
        
        
        selec_all=1;
        for pp=1:size(val_it_all,1)
            pp
            tic
            
            tt=1;
            
            N1_i=size(val_it_all{pp,ind_obj_max},1);
            N2_i=size(val_it_all{pp,interm},1);
            N3_i=size(val_it_all{pp,ind_obj_min},1);
            N1_o=size(val_it_all{pp,ind_obj_max+3},1);
            N2_o=size(val_it_all{pp,interm+3},1);
            N3_o=size(val_it_all{pp,ind_obj_min+3},1);
            ind1_i=[];ind2_i=[];ind3_i=[];
            ind1_o=[];ind2_o=[];ind3_o=[];
            ind1_b=[];ind2_b=[];ind3_b=[];
            for iteri=1:N_shuff
                ind1_i(iteri,:)=randperm(N1_i);
                ind2_i(iteri,:)=randperm(N2_i);
                ind3_i(iteri,:)=randperm(N3_i);
                ind1_o(iteri,:)=randperm(N1_o);
                ind2_o(iteri,:)=randperm(N2_o);
                ind3_o(iteri,:)=randperm(N3_o);
                ind1_b(iteri,:)=randperm(N1_i+N1_o);
                ind2_b(iteri,:)=randperm(N2_i+N2_o);
                ind3_b(iteri,:)=randperm(N3_i+N3_o);
            end
            
            for tt=1:size(time_win,1)
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_both(pp,tt)=mi;
                rate_all_both(pp,1,tt)=nanmean(R1);
                rate_all_both(pp,2,tt)=nanmean(R2);
                rate_all_both(pp,3,tt)=nanmean(R3);
                
                if Shuffle_Preferance
                    R_all=[R1;R2;R3];
                    N1=length(R1);N2=length(R2);N3=length(R3);
                    for iteri=1:N_shuff
                        ind_all=randperm(length(R_all));
                        R1_h=R_all(ind_all(1:N1));
                        R2_h=R_all(ind_all(N1+1:N1+N2));
                        R3_h=R_all(ind_all(end-N3:end));
                        mi=[]; mi= inf_cal(R1_h,R3_h,R2_h);
                        mi_both_shuff(pp,tt,iteri)=mi;
                        rate_both_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_both_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_both_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                    end
                end
                
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_high{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_both(pp,tt)=mi;
                rate_high_both(pp,1,tt)=nanmean(R1);
                rate_high_both(pp,2,tt)=nanmean(R2);
                rate_high_both(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_low{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_both(pp,tt)=mi;
                rate_low_both(pp,1,tt)=nanmean(R1);
                rate_low_both(pp,2,tt)=nanmean(R2);
                rate_low_both(pp,3,tt)=nanmean(R3);
               
                 
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all = [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all = [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all = [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2);nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        ind1=ind1_b(iteri,:);
                        ind2=ind2_b(iteri,:);
                        ind3=ind3_b(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_both_shuff(pp,tt,iteri)=mi;
                        rate_low_both_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_both_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_both_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_both_shuff(pp,tt,iteri)=mi;
                        rate_high_both_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_both_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_both_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_in(pp,tt)=mi;
                rate_all_in(pp,1,tt)=nanmean(R1);
                rate_all_in(pp,2,tt)=nanmean(R2);
                rate_all_in(pp,3,tt)=nanmean(R3);
                
                if Shuffle_Preferance
                    R_all=[R1;R2;R3];
                    N1=length(R1);N2=length(R2);N3=length(R3);
                    for iteri=1:N_shuff
                        ind_all=randperm(length(R_all));
                        R1_h=R_all(ind_all(1:N1));
                        R2_h=R_all(ind_all(N1+1:N1+N2));
                        R3_h=R_all(ind_all(end-N3:end));
                        mi=[]; mi= inf_cal(R1_h,R3_h,R2_h);
                        mi_in_shuff(pp,tt,iteri)=mi;
                        rate_in_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_in_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_in_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                    end
                end
                
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_in(pp,tt)=mi;
                rate_high_in(pp,1,tt)=nanmean(R1);
                rate_high_in(pp,2,tt)=nanmean(R2);
                rate_high_in(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_in(pp,tt)=mi;
                rate_low_in(pp,1,tt)=nanmean(R1);
                rate_low_in(pp,2,tt)=nanmean(R2);
                rate_low_in(pp,3,tt)=nanmean(R3);
                
          
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all= [nanmean(val_it_all{pp,ind_obj_max}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all= [nanmean(val_it_all{pp,interm}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all= [nanmean(val_it_all{pp,ind_obj_min}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        ind1=ind1_i(iteri,:);
                        ind2=ind2_i(iteri,:);
                        ind3=ind3_i(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_in_shuff(pp,tt,iteri)=mi;
                        rate_low_in_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_in_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_in_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_in_shuff(pp,tt,iteri)=mi;
                        rate_high_in_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_in_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_in_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                
                %         end
                %
                %
                %         for tt=1:size(time_win,1)
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_all_out(pp,tt)=mi;
                rate_all_out(pp,1,tt)=nanmean(R1);
                rate_all_out(pp,2,tt)=nanmean(R2);
                rate_all_out(pp,3,tt)=nanmean(R3);
                
                if Shuffle_Preferance
                    R_all=[R1;R2;R3];
                    N1=length(R1);N2=length(R2);N3=length(R3);
                    for iteri=1:N_shuff
                        ind_all=randperm(length(R_all));
                        R1_h=R_all(ind_all(1:N1));
                        R2_h=R_all(ind_all(N1+1:N1+N2));
                        R3_h=R_all(ind_all(end-N3:end));
                        mi=[]; mi= inf_cal(R1_h,R3_h,R2_h);
                        mi_out_shuff(pp,tt,iteri)=mi;
                        rate_out_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_out_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_out_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                    end
                end
                
                
                R1=[]; R2=[]; R3=[];
                R1= [nanmean(val_it_high{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2= [nanmean(val_it_high{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3= [nanmean(val_it_high{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_high_out(pp,tt)=mi;
                rate_high_out(pp,1,tt)=nanmean(R1);
                rate_high_out(pp,2,tt)=nanmean(R2);
                rate_high_out(pp,3,tt)=nanmean(R3);
                
                R1=[]; R2=[]; R3=[];
                R1 = [nanmean(val_it_low{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R2 = [nanmean(val_it_low{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                R3 = [nanmean(val_it_low{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                
                mi=[]; mi= inf_cal(R1,R3,R2);
                mi_low_out(pp,tt)=mi;
                rate_low_out(pp,1,tt)=nanmean(R1);
                rate_low_out(pp,2,tt)=nanmean(R2);
                rate_low_out(pp,3,tt)=nanmean(R3);
                
                if Shuffle_highlow_PPL
                    for iteri=1:N_shuff
                        R1_all=[]; R2_all=[]; R3_all=[];
                        R1_all= [nanmean(val_it_all{pp,ind_obj_max+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R2_all= [nanmean(val_it_all{pp,interm+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        R3_all= [nanmean(val_it_all{pp,ind_obj_min+3}(:,time_win(tt,1):time_win(tt,2)),2)];
                        
                        
                        ind1=ind1_o(iteri,:);
                        ind2=ind2_o(iteri,:);
                        ind3=ind3_o(iteri,:);
                        
                        R1=R1_all(ind1(1:round(Percenti*length(ind1)/100)));
                        R2=R2_all(ind2(1:round(Percenti*length(ind2)/100)));
                        R3=R3_all(ind3(1:round(Percenti*length(ind3)/100)));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_low_out_shuff(pp,tt,iteri)=mi;
                        rate_low_out_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_low_out_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_low_out_shuff(pp,3,tt,iteri)=nanmean(R3);
                        
                        R1=R1_all(ind1(round((100-Percenti)*length(ind1)/100):end));
                        R2=R2_all(ind2(round((100-Percenti)*length(ind2)/100):end));
                        R3=R3_all(ind3(round((100-Percenti)*length(ind3)/100):end));
                        mi=[]; mi= inf_cal(R1,R3,R2);
                        mi_high_out_shuff(pp,tt,iteri)=mi;
                        rate_high_out_shuff(pp,1,tt,iteri)=nanmean(R1);
                        rate_high_out_shuff(pp,2,tt,iteri)=nanmean(R2);
                        rate_high_out_shuff(pp,3,tt,iteri)=nanmean(R3);
                    end
                end
                
            end
            toc
        end
        
%         switch mpi
%             case 1
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
%                     'mi_low_in','mi_low_out','mi_low_both','t_h')
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
%                     'rate_low_in','rate_low_out','rate_low_both','t_h')
%             case 2
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powit_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
%                     'mi_low_in','mi_low_out','mi_low_both','t_h')
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powit_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
%                     'rate_low_in','rate_low_out','rate_low_both','t_h')
%             case 3
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'mi_all_in','mi_all_out','mi_all_both','mi_high_in','mi_high_out','mi_high_both',...
%                     'mi_low_in','mi_low_out','mi_low_both','t_h')
%                 save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'),...
%                     'rate_all_in','rate_all_out','rate_all_both','rate_high_in','rate_high_out','rate_high_both',...
%                     'rate_low_in','rate_low_out','rate_low_both','t_h')
%         end
%         
        if Shuffle_Preferance
             switch mpi
            case 1
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                   'mi_in_shuff','mi_out_shuff','mi_both_shuff')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                     'rate_in_shuff','rate_out_shuff','rate_both_shuff')
            case 2
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_powit_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                   'mi_in_shuff','mi_out_shuff','mi_both_shuff')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_powit_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                     'rate_in_shuff','rate_out_shuff','rate_both_shuff')
            case 3
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                   'mi_in_shuff','mi_out_shuff','mi_both_shuff')
                save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'),...
                     'rate_in_shuff','rate_out_shuff','rate_both_shuff')
        end
           
        end
        
        
        if Shuffle_highlow_PPL
            
            
            save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_wideband_full_shuffled.mat',...
                'mi_high_in_shuff','mi_high_out_shuff','mi_high_both_shuff',...
                'mi_low_in_shuff','mi_low_out_shuff','mi_low_both_shuff',...
                't_h')
        end
        
    end
end
% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_wideband_full_test.mat')
ind_ses=session_select_after_review(1);

%% Load
normalize=0;
mpi=2;
switch mpi
    case 1
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))                  
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_ppl_full_interval_',method,'_',num2str(Percenti),'.mat'))
                   

    case 2
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powit_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powit_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_powit_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_powit_full_interval_',method,'_',num2str(Percenti),'.mat'))
                 
    case 3
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_shuff_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_shuff_powfef_full_interval_',method,'_',num2str(Percenti),'.mat'))
        
end
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


rate_all_in_d= rate_all_in; rate_all_out_d= rate_all_out; rate_all_both_d= rate_all_both;
rate_high_in_d= rate_high_in; rate_high_out_d= rate_high_out; rate_high_both_d= rate_high_both;
rate_low_in_d= rate_low_in; rate_low_out_d= rate_low_out; rate_low_both_d= rate_low_both;

ind_pval=ones(253,1);
load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_interval_300_10ms_wideband_full_shuffled.mat')
if normalize
    %load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_wideband_full.mat')
    
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_wideband_full_shuffled.mat')
    % tri=11
    % %     mi_high_in=squeeze((mi_high_in_shuff(:,:,tri)));mi_low_in=squeeze((mi_low_in_shuff(:,:,tri)));
    % %     mi_high_out=squeeze((mi_high_out_shuff(:,:,tri)));mi_low_out=squeeze((mi_low_out_shuff(:,:,tri)));
    %
    % %      mi_high_in=squeeze(nanmean(mi_high_in_shuff,3));mi_low_in=squeeze(nanmean(mi_low_in_shuff,3));
    % %     mi_high_out=squeeze(nanmean(mi_high_out_shuff,3));mi_low_out=squeeze(nanmean(mi_low_out_shuff,3));
    %
    % %     mi_high_out=mi_high_out_n;mi_low_out=mi_low_out_n;
    % %     mi_high_both=mi_high_both_n;mi_low_both=mi_low_both_n;
    % %
    % %   mi_high_in=(mi_high_in-squeeze(nanmean(mi_high_in_shuff,3)))./squeeze(nanstd(mi_high_in_shuff,[],3));
    % %   mi_low_in=(mi_low_in-squeeze(nanmean(mi_low_in_shuff,3)))./squeeze(nanstd(mi_low_in_shuff,[],3));
    % %   mi_high_out=(mi_high_out-squeeze(nanmean(mi_high_out_shuff,3)))./squeeze(nanstd(mi_high_out_shuff,[],3));
    % %   mi_low_out=(mi_low_out-squeeze(nanmean(mi_low_out_shuff,3)))./squeeze(nanstd(mi_low_out_shuff,[],3));
    %
    %
    % % ind_pval=ones(size(mi_high_in,1),1);
    % %      ind_pval= (nanmean(mi_high_in(:,t_h>700&t_h<950),2)>(squeeze(nanmean(nanmean(mi_high_in_shuff(:,t_h>700&t_h<950,:),2),3))+squeeze(nanstd(nanmean(mi_high_in_shuff(:,t_h>700&t_h<950,:),2),[],3))/sqrt(30)));
    % tt1=900;tt2=950;
    % ind_pval= (abs(nanmean(mi_high_in(:,t_h>tt1&t_h<tt2),2)-...
    %     (squeeze(nanmean(nanmean(mi_high_in_shuff(:,t_h>tt1&t_h<tt2,:),2),3))))>...
    %     (2.045*squeeze(nanstd(nanmean(mi_high_in_shuff(:,t_h>tt1&t_h<tt2,:),2),[],3))/sqrt(size(mi_high_in_shuff,3))));
    %
    % %      ind_pval=ind_pval& (nanmean(mi_low_in(:,t_h>700&t_h<950),2)>((squeeze(nanmean(nanmean(mi_low_in_shuff(:,t_h>700&t_h<950,:),2),3))+squeeze(nanstd(nanmean(mi_low_in_shuff(:,t_h>700&t_h<950,:),2),[],3))/sqrt(30))));
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_interval_300_10ms_wideband_full.mat')
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_interval_300_10ms_wideband_full_shuffled.mat')
    %
    
    %
    % ind_pval= (abs(mi_high_in-...
    %     (squeeze(nanmean(mi_high_in_shuff,3))))>...
    %     (1.962*squeeze(nanstd(mi_high_in_shuff,[],3)))/sqrt(size(mi_high_in_shuff,3)));
    %
    ind_pval1= mi_high_in  >  squeeze(prctile(mi_high_in_shuff,80,3));
    ind_pval2= mi_high_in  <  squeeze(prctile(mi_high_in_shuff,20,3));
    ind_pval=ind_pval2;
end

switch mpi
    case 1
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_ppl_full_',method,'_',num2str(Percenti),'.mat'))
    case 2
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powit_full_',method,'_',num2str(Percenti),'.mat'))
    case 3
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\ROC_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat'))
        load(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\Rate_timeCoures_300_10ms_powfef_full_',method,'_',num2str(Percenti),'.mat'))
end
load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\rate It\plv_categorized_it_rate_wavelet_wideband_full.mat')

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

ix_h_tri_in=(trail_num_in>3)';
ix_h_tri_out=(trail_num_out>3)';

% ix_h_tri_in=(trail_num_in>3)'&in_rf';
% ix_h_tri_out=(trail_num_out>3)'&out_rf';

  ix_h_tri_in=ones(253,1); ix_h_tri_out=ones(253,1);

%% Session selection
interval=1;
M1=0;
M2=0;
tresh=-1;
status = 1;   % status = 1 All , status = 2 SFN; status = 3 performance
ind_ses=session_select_after_review(status);
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
pp=0;ind_it_neurons_ses=[];ind_it_p=[];pval_dsi_it=[];
for ss= 1:92
    lfp=gen_lfp_pairs(ss);
    
    
    
    
    IT_=resp(ss).IT;
    if ~isempty(IT_)
        for nn=1:size(IT_,2)
            IT=IT_{nn};
            pp=pp+1;
            ind_it_neurons_ses=[ind_it_neurons_ses ss];
            
            obj_ix_h=[resp(ss).preff_obj(nn) resp(ss).npreff_obj(nn)];
            obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
            
            
            
            
            ix_h = [];
            ix_h = ismember(resp(ss).condition{:},obj_ix_h(1)+[0 3])&(resp(ss).bhv{:}==1);
            val_1 = [];  val_1 = signrank(mean(IT(ix_h,1100:1700),2),mean(IT(ix_h,1:500),2));
            
            ix_h = [];
            ix_h = ismember(resp(ss).condition{:},obj_ix_h(2)+[0 3])&(resp(ss).bhv{:}==1);
            val_2 = [];  val_2 = signrank(mean(IT(ix_h,1100:1700),2),mean(IT(ix_h,1:500),2));
            
            ix_h = [];
            ix_h = ismember(resp(ss).condition{:},obj_ix_h(3)+[0 3])&(resp(ss).bhv{:}==1);
            val_3 = [];  val_3 = signrank(mean(IT(ix_h,1100:1700),2),mean(IT(ix_h,1:500),2));
            
            ind_it_p=[ind_it_p (val_1<0.05)|(val_2<0.05)|(val_3<0.05)];
            
            %         sig=[];  sig_=[];
            %         sig=resp(ss_p).IT{nn}(:,1:2100);
            %         N_trials=size(sig,1);
            %         for tt=1:size(sig,1)
            %             sig_(tt,:)=SmoothN(sig(tt,:)*1, 10);
            %         end
            %         N_bts=1001;
            %         fr_in= mean2(sig_(logical(ismember(resp(ss_p).condition{:},[obj_ix_h(1) obj_ix_h(1)+3])),t1:t2));
            %         fr_out= mean2(sig_(logical(ismember(resp(ss_p).condition{:},[obj_ix_h(3) obj_ix_h(3)+3])),t1:t2));
            %         si=(fr_in-fr_out)/(fr_in+fr_out);
            %         for iteri=1:N_bts
            %             ind_h=randperm(N_trials,round(N_trials/3));
            %             fr_in_h= mean2(sig_(ind_h,t1:t2));
            %
            %             ind_h_2=setdiff([1:N_trials],ind_h);
            %             ind_h_ind=randperm(length(ind_h_2),round(N_trials/3));
            %             ind_h=ind_h_2(ind_h_ind);
            %             fr_out_h= mean2(sig_(ind_h,t1:t2));
            %             si_h(iteri)=(fr_in_h-fr_out_h)/(fr_in_h+fr_out_h);
            %         end
            %         pval_dsi_it(pp)=1-sum(si>si_h)/N_bts;
            %
            
            
        end
    end
end
ind_it_neurons=ismember(ind_it_neurons_ses,ind_ses);

%%Neurons selection
close all
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

figure('name','ROC IN')
subplot(121)
hold on
x_h_=[];x_h_=mi_low_n_(ix_h_in);
y_h_=[];y_h_=mi_high_n_(ix_h_in);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);

eval(strcat('selectivity_low_high_in= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

scatter(x_h,y_h,'k')
text( 0.5,0.5,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.5,0.4,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')

set(gca,'fontsize',14,'fontweight','bold')
xlabel('AUC Low PPL')
ylabel('AUC High PPL')
 axis([0.2 1.1 0.2 1.1])
% axis([-0.5 1.1 -0.5 1.1])

%  xbin=[-0.5:0.03:0.5];
subplot(122)
var_=[]; var_=-1*(y_h-x_h);
% hist(var_,bin_num)
bin_=0.05;
nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
[counts,centers]=hist(var_,nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)

m=median(var_);
line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
line([m m],ylim,'linewidth',2,'color','r')
if(p<0.05)
    text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
    
    text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
end
set(gca,'fontsize',14,'fontweight','bold')
ylabel('#')
xlim([-0.4 0.4])

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

figure('name','OUT')
subplot(121)
hold on
x_h_=[];x_h_=mi_low_n_(ix_h_out);
y_h_=[];y_h_=mi_high_n_(ix_h_out);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
plot([-0.5:.1:1.5],[-0.5:0.1:1.5])
p=signrank(x_h,y_h)
n=length(y_h);
eval(strcat('selectivity_low_high_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))
scatter(x_h,y_h,'k')
text( 0.5,0.5,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.5,0.4,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
xlabel('AUC Low PPL')
ylabel('AUC High PPL')
 axis([0.2 1.1 0.2 1.1])
% axis([-0.5 1.1 -0.5 1.1])


xbin=[-0.5:0.03:0.5];
subplot(122)
var_=[]; var_=-1*(y_h-x_h);
bin_=0.05;
nbin1=-2:bin_:2.5;nbin2=-0.5:bin_:0.5;
[counts,centers]=hist(var_,nbin2);
counts=counts/sum(counts);
bar(centers,counts,1)
line([0 0],ylim)

m=median(var_);
line([0 0],ylim,'linestyle',':','linewidth',2,'color','k')
line([m m],ylim,'linewidth',2,'color','r')
if(p<0.05)
    text(m,5,num2str(strcat('med = ',num2str(m,2),';')),'fontsize',14,'fontweight','bold')
    
    text(m,4,num2str(strcat('p = ',num2str(p,2),';')),'fontsize',14,'fontweight','bold')
end
set(gca,'fontsize',14,'fontweight','bold')
ylabel('#')
xlim([-0.4 0.4])

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preveiw time course roc
bb_y_low=0.47;bb_y_high=0.85;

win=1;
figure ('name','IN')
ax_=subplot(1,2,1);
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



win=1;

ax_=subplot(1,2,2);
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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MI time coures
bb_y_low=-0.1;bb_y_high=0.1;


win=1;
figure ('name','IN')
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
figure
subplot(121)

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

eval(strcat('selectivity_mi_in_out= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')


text( 0.2,0.3,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.2,0.2,strcat('n = ',num2str(size(x_h,1),3)),'fontsize',14,'fontweight','bold')

set(gca,'fontsize',14,'fontweight','bold')
xlabel('MI OUT')
ylabel('MI IN')
 axis([-0.4 0.6 -0.4 0.6])
subplot(122)
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

%%Fig behrad
if ~interval
rate_high_in_norm=rate_high_in;rate_high_out_norm=rate_high_out;
rate_low_in_norm=rate_low_in;rate_low_out_norm=rate_low_out;
normalized=0;
a=1;
if normalized
    a=1;
    
    for pp=1:size(rate_high_in,1)
        sig_=[];  sig_=[squeeze(rate_high_in(pp,:,:));squeeze(rate_high_out(pp,:,:));...
            squeeze(rate_low_in(pp,:,:));squeeze(rate_low_out(pp,:,:))];
        var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        var_max=[];  var_max=max(max(sig_));
        
        %  sig_=[];  sig_=[squeeze(rate_high_in(pp,:,:));...
        %         squeeze(rate_low_in(pp,:,:))];
        %       var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=[];  var_max=max(max(sig_));
        
        %            sig_=[];
        %              sig_=[squeeze(rate_high_in(pp,:,:))];
        %      var_b=[]; var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %           var_max=[];   var_max=max(max(sig_));
        rate_high_in_norm(pp,:,:)=(rate_high_in(pp,:,:)-var_b)./(var_max-var_b);
        
        %
        %                 sig_=[];  sig_=[squeeze(rate_low_in(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_low_in_norm(pp,:,:)=(rate_low_in(pp,:,:)-var_b)./(var_max-var_b);
        
        %           sig_=[];  sig_=[squeeze(rate_high_out(pp,:,:));...
        %         squeeze(rate_low_out(pp,:,:))];
        %       var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=[];  var_max=max(max(sig_));
        %
        
        %                          sig_=[];     sig_=[squeeze(rate_high_out(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_high_out_norm(pp,:,:)=(rate_high_out(pp,:,:)-var_b)./(var_max-var_b);
        
        %
        %                 sig_=[];  sig_=[squeeze(rate_low_out(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_low_out_norm(pp,:,:)=(rate_low_out(pp,:,:)-var_b)./(var_max-var_b);
        
    end
    
end
%%% IN
inPrefHigh=a*squeeze(nanmean(rate_high_in_norm(ix_h_in,1,t_h>t1_delay&t_h<t2_delay),3));
inNPrefHigh=a*squeeze(nanmean(rate_high_in_norm(ix_h_in,3,t_h>t1_delay&t_h<t2_delay),3));
inPrefLow=a*squeeze(nanmean(rate_low_in_norm(ix_h_in,1,t_h>t1_delay&t_h<t2_delay),3));
inNPrefLow=a*squeeze(nanmean(rate_low_in_norm(ix_h_in,3,t_h>t1_delay&t_h<t2_delay),3));
outPrefHigh=a*squeeze(nanmean(rate_high_out_norm(ix_h_out,1,t_h>t1_delay&t_h<t2_delay),3));
outNPrefHigh=a*squeeze(nanmean(rate_high_out_norm(ix_h_out,3,t_h>t1_delay&t_h<t2_delay),3));
outPrefLow=a*squeeze(nanmean(rate_low_out_norm(ix_h_out,1,t_h>t1_delay&t_h<t2_delay),3));
outNPrefLow=a*squeeze(nanmean(rate_low_out_norm(ix_h_out,3,t_h>t1_delay&t_h<t2_delay),3));
else
 
rate_high_in_norm=rate_high_in_d;rate_high_out_norm=rate_high_out_d;
rate_low_in_norm=rate_low_in_d;rate_low_out_norm=rate_low_out_d;
normalized=0;
a=1;
%%% IN
inPrefHigh=a*rate_high_in_norm(ix_h_in,1);
inNPrefHigh=a*rate_high_in_norm(ix_h_in,3);
inPrefLow=a*rate_low_in_norm(ix_h_in,1);
inNPrefLow=a*rate_low_in_norm(ix_h_in,3);
outPrefHigh=a*rate_high_out_norm(ix_h_out,1);
outNPrefHigh=a*rate_high_out_norm(ix_h_out,3);
outPrefLow=a*rate_low_out_norm(ix_h_out,1);
outNPrefLow=a*rate_low_out_norm(ix_h_out,3);   
end

sig_=[];sig_=[inPrefHigh inNPrefHigh inPrefLow inNPrefLow];
max_=repmat(max(sig_,[],2),1,4);base_=repmat(min(sig_,[],2),1,4);
sig=(sig_-base_)./(max_-base_);
inPrefHigh=sig(:,1);inNPrefHigh=sig(:,2);inPrefLow=sig(:,3);inNPrefLow=sig(:,4);
figure('name','selected')
ax=subplot(121);
hold on
plot(nanmean([inNPrefHigh inPrefHigh]),'g')
plot(nanmean([inNPrefLow inPrefLow]),'k')

errorbar(nanmean([inNPrefHigh inPrefHigh]),nanstd([inNPrefHigh inPrefHigh])./sqrt(length(inNPrefHigh)),'g')
errorbar(nanmean([inNPrefLow inPrefLow]),nanstd([inNPrefLow inPrefLow])./sqrt(length(inNPrefLow)),'k')
ax.XTick=[1 2];
ax.XTickLabel=[{'npref'} {'pref'}];

ylabel('Firing rate ')
ylim([0 1])
xlim([0.8 2.1])
set(gca,'fontsize',14,'fontweight','bold')
title('IN')
y_h=inPrefHigh;
x_h=inPrefLow;
n=size(inPrefHigh,1);
p=signrank(x_h,y_h)
eval(strcat('fr_low_high_in_preff= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

y_h=inNPrefHigh;
x_h=inNPrefLow;
n=size(inNPrefHigh,1);
p=signrank(x_h,y_h)
eval(strcat('fr_low_high_in_npreff= round([nanmean(y_h) , nanmean(x_h), nanmean(y_h)-nanmean(x_h),nanstd(y_h-x_h)/sqrt(n),n, p],3)'))

%%% OUT

sig_=[];sig_=[outPrefHigh outNPrefHigh outPrefLow outNPrefLow];
max_=repmat(max(sig_,[],2),1,4);base_=repmat(min(sig_,[],2),1,4);
sig=(sig_-base_)./(max_-base_);
outPrefHigh=sig(:,1);outNPrefHigh=sig(:,2);outPrefLow=sig(:,3);outNPrefLow=sig(:,4);

ax=subplot(122);
hold on
plot(nanmean([outNPrefHigh outPrefHigh]),'g')
plot(nanmean([outNPrefLow outPrefLow]),'k')

errorbar(nanmean([outNPrefHigh outPrefHigh]),nanstd([outNPrefHigh outPrefHigh])./sqrt(length(outNPrefHigh)),'g')
errorbar(nanmean([outNPrefLow outPrefLow]),nanstd([outNPrefLow outPrefLow])./sqrt(length(outNPrefLow)),'k')
ylabel('Firing rate ')
ylim([0 1])
xlim([0.8 2.1])
ax.XTick=[1 2];
ax.XTickLabel=[{'npref'} {'pref'}];
set(gca,'fontsize',14,'fontweight','bold')
title('OUT')

%%Fig behrad for all
if ~interval
rate_high_in_norm=rate_high_in;rate_high_out_norm=rate_high_out;
rate_low_in_norm=rate_low_in;rate_low_out_norm=rate_low_out;
normalized=0;
a=1;
if normalized
    a=1;
    
    for pp=1:size(rate_high_in,1)
        sig_=[];  sig_=[squeeze(rate_high_in(pp,:,:));squeeze(rate_high_out(pp,:,:));...
            squeeze(rate_low_in(pp,:,:));squeeze(rate_low_out(pp,:,:))];
        var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        var_max=[];  var_max=max(max(sig_));
        
        %  sig_=[];  sig_=[squeeze(rate_high_in(pp,:,:));...
        %         squeeze(rate_low_in(pp,:,:))];
        %       var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=[];  var_max=max(max(sig_));
        
        %            sig_=[];
        %              sig_=[squeeze(rate_high_in(pp,:,:))];
        %      var_b=[]; var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %           var_max=[];   var_max=max(max(sig_));
        rate_high_in_norm(pp,:,:)=(rate_high_in(pp,:,:)-var_b)./(var_max-var_b);
        
        %
        %                 sig_=[];  sig_=[squeeze(rate_low_in(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_low_in_norm(pp,:,:)=(rate_low_in(pp,:,:)-var_b)./(var_max-var_b);
        
        %           sig_=[];  sig_=[squeeze(rate_high_out(pp,:,:));...
        %         squeeze(rate_low_out(pp,:,:))];
        %       var_b=[];  var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=[];  var_max=max(max(sig_));
        %
        
        %                          sig_=[];     sig_=[squeeze(rate_high_out(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_high_out_norm(pp,:,:)=(rate_high_out(pp,:,:)-var_b)./(var_max-var_b);
        
        %
        %                 sig_=[];  sig_=[squeeze(rate_low_out(pp,:,:))];
        %       var_b=(mean2(sig_(:,t_h>-500&t_h<0)));
        %             var_max=max(max(sig_));
        rate_low_out_norm(pp,:,:)=(rate_low_out(pp,:,:)-var_b)./(var_max-var_b);
        
    end
    
end
%%% IN
inPrefHigh=a*squeeze(nanmean(rate_high_in_norm(:,1,t_h>t1_delay&t_h<t2_delay),3));
inNPrefHigh=a*squeeze(nanmean(rate_high_in_norm(:,3,t_h>t1_delay&t_h<t2_delay),3));
inPrefLow=a*squeeze(nanmean(rate_low_in_norm(:,1,t_h>t1_delay&t_h<t2_delay),3));
inNPrefLow=a*squeeze(nanmean(rate_low_in_norm(:,3,t_h>t1_delay&t_h<t2_delay),3));
outPrefHigh=a*squeeze(nanmean(rate_high_out_norm(:,1,t_h>t1_delay&t_h<t2_delay),3));
outNPrefHigh=a*squeeze(nanmean(rate_high_out_norm(:,3,t_h>t1_delay&t_h<t2_delay),3));
outPrefLow=a*squeeze(nanmean(rate_low_out_norm(:,1,t_h>t1_delay&t_h<t2_delay),3));
outNPrefLow=a*squeeze(nanmean(rate_low_out_norm(:,3,t_h>t1_delay&t_h<t2_delay),3));
else
 
rate_high_in_norm=rate_high_in_d;rate_high_out_norm=rate_high_out_d;
rate_low_in_norm=rate_low_in_d;rate_low_out_norm=rate_low_out_d;
normalized=0;
a=1;
%%% IN
inPrefHigh=a*rate_high_in_norm(:,1);
inNPrefHigh=a*rate_high_in_norm(:,3);
inPrefLow=a*rate_low_in_norm(:,1);
inNPrefLow=a*rate_low_in_norm(:,3);
outPrefHigh=a*rate_high_out_norm(:,1);
outNPrefHigh=a*rate_high_out_norm(:,3);
outPrefLow=a*rate_low_out_norm(:,1);
outNPrefLow=a*rate_low_out_norm(:,3);   
end

sig_=[];sig_=[inPrefHigh inNPrefHigh inPrefLow inNPrefLow outPrefHigh outNPrefHigh outPrefLow outNPrefLow];
max_=repmat(max(sig_,[],2),1,8);base_=repmat(min(sig_,[],2),1,8);
sig=(sig_-base_)./(max_-base_);
inPrefHigh=sig(:,1);inNPrefHigh=sig(:,2);inPrefLow=sig(:,3);inNPrefLow=sig(:,4);
outPrefHigh=sig(:,5);outNPrefHigh=sig(:,6);outPrefLow=sig(:,7);outNPrefLow=sig(:,8);

figure('name','All')
ax=subplot(121);
hold on
plot(nanmean([inNPrefHigh inPrefHigh]),'g')
plot(nanmean([inNPrefLow inPrefLow]),'k')

errorbar(nanmean([inNPrefHigh inPrefHigh]),nanstd([inNPrefHigh inPrefHigh])./sqrt(length(inNPrefHigh)),'g')
errorbar(nanmean([inNPrefLow inPrefLow]),nanstd([inNPrefLow inPrefLow])./sqrt(length(inNPrefLow)),'k')
ax.XTick=[1 2];
ax.XTickLabel=[{'npref'} {'pref'}];

ylabel('Firing rate ')
%  ylim([-0.2 0.1])
xlim([0.8 2.1])
set(gca,'fontsize',14,'fontweight','bold')
title('IN')
%%% OUT
% outPrefHigh=a*squeeze(nanmean(rate_high_out_norm(ix_h_out,1,t_h>t1&t_h<t2),3));
% outNPrefHigh=a*squeeze(nanmean(rate_high_out_norm(ix_h_out,3,t_h>t1&t_h<t2),3));
% outPrefLow=a*squeeze(nanmean(rate_low_out_norm(ix_h_out,1,t_h>t1&t_h<t2),3));
% outNPrefLow=a*squeeze(nanmean(rate_low_out_norm(ix_h_out,3,t_h>t1&t_h<t2),3));
ax=subplot(122);
hold on
plot(nanmean([outNPrefHigh outPrefHigh]),'g')
plot(nanmean([outNPrefLow outPrefLow]),'k')

errorbar(nanmean([outNPrefHigh outPrefHigh]),nanstd([outNPrefHigh outPrefHigh])./sqrt(length(outNPrefHigh)),'g')
errorbar(nanmean([outNPrefLow outPrefLow]),nanstd([outNPrefLow outPrefLow])./sqrt(length(outNPrefLow)),'k')
ylabel('Firing rate ')
%  ylim([-0.2 0.1])
xlim([0.8 2.1])
ax.XTick=[1 2];
ax.XTickLabel=[{'npref'} {'pref'}];
set(gca,'fontsize',14,'fontweight','bold')
title('OUT')

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

figure('name','IN')
subplot(121)
title('IN')
hold on
% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')

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
subplot(122)
% figure('name','OUT')
title('OUT')
hold on
% scatter(x_h,y_h)
scatter(x_h(m1_ind),y_h(m1_ind),'b*')
scatter(x_h(m2_ind),y_h(m2_ind),'gs')

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

end