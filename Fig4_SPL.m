function Fig4_SPL(resp)

%% Set Paths
% clear all
% clc
% close all

%% Settings
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
M1=0;
M2=0;

normalize_= @(x) x/nansum(x);
plv=@(x) abs(nanmean(exp(1j*x)));
control=0;
% freqs2use  =[4:1:70 71:5:130];
 freqs2use  = 4:2:60;
%freqs2use= [3 4 10 12 15];

Color_loc=[1,0,0;0,0,1];
Color_obj=[1,0.5,0;0.5,0.5,0.5;0,1,0.5];
Color_bhv=[1,1,0;0,1,1];
win_=1;
t1=1100;t2=1700;N_bts=1001;

noramlized=0;


% freq_bands=[3.9,7.9,19.9,34.9,49.9;8.1,15.1,28.1,50.1,130.1];

freq_bands=[3.9,7.9,21.9,34.9,49.9,49.9;8.1,15.1,28.1,50.1,130.1,150.1];

%%
gen_lfp_pairs(1).IT=resp(1).LFP_IT;
gen_lfp_pairs(1).FEF=resp(1).LFP_FEF;
gen_lfp_pairs(1).bhv{1}=resp(1).bhv;gen_lfp_pairs(1).condition{1}=resp(1).condition';
gen_lfp_pairs(1).preff_obj=resp(1).pref_obj_mu;
gen_lfp_pairs(1).npreff_obj=resp(1).npref_obj_mu;
gen_lfp_pairs(1).FEF_artifact=resp(1).FEF_artifact;
gen_lfp_pairs(1).IT_artifact=resp(1).IT_artifact;


%% Across area calculation 
%%Calculate
disp('Wait! Calculating SPL  **********')

%%%  parfor Methode 2 plv in each condition trails matched
cal=1;
if cal
    
    st_ind=1;end_ind=2100;
    
    plv_cr=[];phase_fef_cr=[];phase_it_cr=[];
    plv_wr=[];phase_fef_wr=[];phase_it_wr=[];
    t1_b=-500;t2_b=0;
    t1_ld=600;t2_ld=1200;
    
    min_RT=50;
    n_hist_bins=40;
    xbins =linspace(-pi,pi,n_hist_bins+1);
    nn=0;nn=0;it_tag=[];fef_tag=[];
    N_shuff=30;N_sub=100;
    Info_aa=[];
    for ss=1:size(gen_lfp_pairs,2)
        
        
        tic
        lfp=gen_lfp_pairs(ss);
        trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
        
        s_it = []; s_it = lfp.IT(trii_h,st_ind:end_ind);
        s_fef = []; s_fef = lfp.FEF(trii_h,st_ind:end_ind);
        
        t=1:2100;
        sig_it = [];sig_fef = [];
        [sig_it,sig_fef,freqs2use_,t]= wavelet_core_pertrial(s_it,s_fef,freqs2use);
        %  [sig_it,sig_fef]= hilbert_core_map(s_it,s_fef);
        t_h=t-500;
        phase_fef_ = [];phase_fef_ = angle(sig_fef);
        phase_it_ = [];phase_it_ = angle(sig_it);
        diff_phi_ =[];diff_phi_=angle(exp(1j*(angle(sig_fef)-angle(sig_it))));
        trial_num=size(diff_phi_,1);
        
        val_shuff_it_cr = [];
        val_shuff_it_wr = [];
        val_shuff_fef_cr = [];
        val_shuff_fef_wr = [];
        val_fef_it_cr = [];
        val_it_it_cr = [];
        val_diff_it_cr = [];
        val_fef_it_n_cr = [];
        val_it_it_n_cr = [];
        val_diff_it_n_cr = [];
        val_fef_it_wr = [];
        val_it_it_wr = [];
        val_diff_it_wr = [];
        val_fef_it_n_wr = [];
        val_it_it_n_wr = [];
        val_diff_it_n_wr = [];
        N_sam_sp=30;
        t_h_sp=t-500;
        IT=[];
        
        IT=resp(ss).IT;
        if ss==92
            IT=resp(ss).FEF{1};
            resp(ss).preff_obj(1)=1;
            resp(ss).npreff_obj(1)=3;
        end
        
        for nn=1:size(IT,2)
            
            obj_ix_h=[resp(ss).pref_obj_su(nn) resp(ss).npref_obj_su(nn)];
            obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
            it_tag=[it_tag ss*10+nn];
            it_=IT{nn}(trii_h,:);
            
            rate_it= (it_(:,t1_ld<t_h_sp&t_h_sp<t2_ld));
            n_trials=size(rate_it,1);
            for fi=1:length(freqs2use)
                
                val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
                phase_fef =squeeze(phase_fef_(:,fi,t1_ld<t_h&t_h<t2_ld));
                phase_it =squeeze(phase_it_(:,fi,t1_ld<t_h&t_h<t2_ld));
                diff_phi=squeeze(diff_phi_(:,fi,t1_ld<t_h&t_h<t2_ld));
                for condi=1:6
                    if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
                    ix_h_cr = [];
                    ix_h_cr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==1)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    ix_h_wr = [];
                    ix_h_wr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==0)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    
                    N_cr=sum(ix_h_cr);N_wr=sum(ix_h_wr);
                    
                    phase_fef_cr = []; phase_fef_cr=phase_fef(ix_h_cr,:);
                    phase_it_cr = [];phase_it_cr=phase_it(ix_h_cr,:);
                    diff_phi_cr = [];diff_phi_cr=diff_phi(ix_h_cr,:);
                    rate_it_cr = [];rate_it_cr = rate_it(ix_h_cr,:);
                    val_fef=[];val_it=[]; val_shuff=[];val_shuff=[];
                    [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_cr,phase_it_cr,rate_it_cr,N_sam_sp);
                    
                    %      [coho,scoho,phaso,urcoho,srcoho,ifcoher] = compute_fries_coherence_match(s_it(ix_h_cr,t1_ld<t_h&t_h<t2_ld),rate_it_cr,100,N_sam_sp)
                    
                    val_fef_it_cr(nn,condi,fi)=(val_fef);
                    val_it_it_cr(nn,condi,fi)=(val_it);
                    val_shuff_it_cr(nn,condi,fi)=(val_shuff);
                    
                    if N_wr>0
                        phase_fef_wr = []; phase_fef_wr=phase_fef(ix_h_wr,:);
                        phase_it_wr = [];phase_it_wr=phase_it(ix_h_wr,:);
                        diff_phi_wr = [];diff_phi_wr=diff_phi(ix_h_wr,:);
                        rate_it_wr = [];rate_it_wr = rate_it(ix_h_wr,:);
                        val_fef=[];val_it=[]; val_shuff=[];
                        [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_wr,phase_it_wr,rate_it_wr,N_sam_sp);
                        
                        val_fef_it_wr(nn,condi,fi)=(val_fef);
                        val_it_it_wr(nn,condi,fi)=(val_it);
                        val_shuff_it_wr(nn,condi,fi)=(val_shuff);
                    else
                        val_fef_it_wr(nn,condi,fi)=nan;
                        val_it_it_wr(nn,condi,fi)=nan;
                        val_shuff_it_wr(nn,condi,fi)=nan;
                    end
                    
                    
                end
            end
        end
        Info_aa(ss).val_fef_it_cr=val_fef_it_cr;
        Info_aa(ss).val_it_it_cr=val_it_it_cr;
        Info_aa(ss).val_shuff_it_cr=val_shuff_it_cr;
        
        Info_aa(ss).val_fef_it_wr=val_fef_it_wr;
        Info_aa(ss).val_it_it_wr=val_it_it_wr;
        Info_aa(ss).val_shuff_it_wr=val_shuff_it_wr;
        
        val_shuff_fef_cr = [];
        val_shuff_fef_wr = [];
        
        val_fef_fef_cr = [];
        val_it_fef_cr = [];
        val_diff_fef_cr = [];
        val_fef_fef_n_cr = [];
        val_it_fef_n_cr = [];
        val_diff_fef_n_cr = [];
        
        val_fef_fef_wr = [];
        val_it_fef_wr = [];
        val_diff_fef_wr = [];
        val_fef_fef_n_wr = [];
        val_it_fef_n_wr = [];
        val_diff_fef_n_wr = [];
        obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        FEF=resp(ss).FEF;
        for nn=1:size(FEF,2)
            
            fef_tag=[fef_tag ss*10+nn];
            fef_=FEF{nn}(trii_h,:);
            
            rate_fef= (fef_(:,t1_ld<t_h_sp&t_h_sp<t2_ld));
            n_trials=size(rate_fef,1);
            for fi=1:length(freqs2use)
                
                val_fef_fef_=[];val_it_fef_=[];val_diff_fef_=[];
                phase_fef =squeeze(phase_fef_(:,fi,t1_ld<t_h&t_h<t2_ld));
                phase_it =squeeze(phase_it_(:,fi,t1_ld<t_h&t_h<t2_ld));
                diff_phi=squeeze(diff_phi_(:,fi,t1_ld<t_h&t_h<t2_ld));
                for condi=1:6
                    if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
                    ix_h_cr = [];
                    ix_h_cr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==1)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    ix_h_wr = [];
                    ix_h_wr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==0)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    
                    N_cr=sum(ix_h_cr);N_wr=sum(ix_h_wr);
                    
                    phase_fef_cr = []; phase_fef_cr=phase_fef(ix_h_cr,:);
                    phase_it_cr = [];phase_it_cr=phase_it(ix_h_cr,:);
                    diff_phi_cr = [];diff_phi_cr=diff_phi(ix_h_cr,:);
                    rate_fef_cr = [];rate_fef_cr = rate_fef(ix_h_cr,:);
                    val_fef=[];val_it=[];val_diff=[]; val_fef_n=[];val_it_n=[];val_diff_n=[];
                    [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_cr,phase_it_cr,rate_fef_cr,N_sam_sp);
                    
                    val_fef_fef_cr(nn,condi,fi)=(val_fef);
                    val_it_fef_cr(nn,condi,fi)=(val_it);
                    val_shuff_fef_cr(nn,condi,fi)=(val_shuff);
                    
                    if N_wr>0
                        phase_fef_wr = []; phase_fef_wr=phase_fef(ix_h_wr,:);
                        phase_it_wr = [];phase_it_wr=phase_it(ix_h_wr,:);
                        diff_phi_wr = [];diff_phi_wr=diff_phi(ix_h_wr,:);
                        rate_fef_wr = [];rate_fef_wr = rate_fef(ix_h_wr,:);
                        val_fef=[];val_it=[];val_diff=[]; val_fef_n=[];val_it_n=[];val_diff_n=[];
                        [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_wr,phase_it_wr,rate_fef_wr,N_sam_sp);
                        
                        val_fef_fef_wr(nn,condi,fi)=(val_fef);
                        val_it_fef_wr(nn,condi,fi)=(val_it);
                        val_shuff_fef_wr(nn,condi,fi)=(val_shuff);
                        
                    else
                        val_fef_fef_wr(nn,condi,fi)=nan;
                        val_it_fef_wr(nn,condi,fi)=nan;
                        val_shuff_fef_wr(nn,condi,fi)=nan;
                        
                        
                    end
                    
                    
                    
                    
                    
                end
                
            end
        end
        
        Info_aa(ss).val_fef_fef_cr=val_fef_fef_cr;
        Info_aa(ss).val_it_fef_cr=val_it_fef_cr;
        Info_aa(ss).val_shuff_fef_cr=val_shuff_fef_cr;
        
        Info_aa(ss).val_fef_fef_wr=val_fef_fef_wr;
        Info_aa(ss).val_it_fef_wr=val_it_fef_wr;
        Info_aa(ss).val_shuff_fef_wr=val_shuff_fef_wr;
        
        ss
        toc
    end
    
    
%     save('D:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\spl_ppc_per_session_bahmani_full.mat'...
%         , 'Info','freqs2use')
end
%% Within area calculation 
%%Calculate

%%%  parfor Methode 2 plv in each condition trails matched
cal=1;
if cal
    
    st_ind=1;end_ind=2100;
    
    plv_cr=[];phase_fef_cr=[];phase_it_cr=[];
    plv_wr=[];phase_fef_wr=[];phase_it_wr=[];
    t1_b=-500;t2_b=0;
    t1_ld=600;t2_ld=1200;
    
    min_RT=50;
    n_hist_bins=40;
    xbins =linspace(-pi,pi,n_hist_bins+1);
    nn=0;nn=0;it_tag=[];fef_tag=[];
    N_shuff=30;N_sub=100;
    Info_wa=[];
    for ss=1:size(gen_lfp_pairs,2)
        
        
        tic
        lfp=gen_lfp_pairs(ss);
        trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
        
        s_it = []; s_it = lfp.IT(trii_h,st_ind:end_ind);
        s_fef = []; s_fef = lfp.FEF(trii_h,st_ind:end_ind);
        
        t=1:2100;
        sig_it = [];sig_fef = [];
        [sig_it,sig_fef,freqs2use_,t]= wavelet_core_pertrial(s_it,s_fef,freqs2use);
        %  [sig_it,sig_fef]= hilbert_core_map(s_it,s_fef);
        t_h=t-500;
        phase_fef_ = [];phase_fef_ = angle(sig_fef);
        phase_it_ = [];phase_it_ = angle(sig_it);
        diff_phi_ =[];diff_phi_=angle(exp(1j*(angle(sig_fef)-angle(sig_it))));
        trial_num=size(diff_phi_,1);
        
        val_shuff_it_cr = [];
        val_shuff_it_wr = [];
        val_shuff_fef_cr = [];
        val_shuff_fef_wr = [];
        val_fef_it_cr = [];
        val_it_it_cr = [];
        val_diff_it_cr = [];
        val_fef_it_n_cr = [];
        val_it_it_n_cr = [];
        val_diff_it_n_cr = [];
        val_fef_it_wr = [];
        val_it_it_wr = [];
        val_diff_it_wr = [];
        val_fef_it_n_wr = [];
        val_it_it_n_wr = [];
        val_diff_it_n_wr = [];
        N_sam_sp=30;
        t_h_sp=t-500;
        IT=[];
        
        IT=resp(ss).IT;
        if ss==92
            IT=resp(ss).FEF{1};
            resp(ss).preff_obj(1)=1;
            resp(ss).npreff_obj(1)=3;
        end
        
        for nn=1:size(IT,2)
            
            obj_ix_h=[resp(ss).pref_obj_su(nn) resp(ss).npref_obj_su(nn)];
            obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
            it_tag=[it_tag ss*10+nn];
            it_=IT{nn}(trii_h,:);
            
            rate_it= (it_(:,t1_ld<t_h_sp&t_h_sp<t2_ld));
            n_trials=size(rate_it,1);
            for fi=1:length(freqs2use)
                
                val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
                phase_fef =squeeze(phase_fef_(:,fi,t1_ld<t_h&t_h<t2_ld));
                phase_it =squeeze(phase_it_(:,fi,t1_ld<t_h&t_h<t2_ld));
                diff_phi=squeeze(diff_phi_(:,fi,t1_ld<t_h&t_h<t2_ld));
                for condi=1:3
                    condi_h=[obj_ix_h(condi) obj_ix_h(condi)+3];
                    ix_h_cr = [];
                    ix_h_cr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==1)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    ix_h_wr = [];
                    ix_h_wr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==0)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    
                    N_cr=sum(ix_h_cr);N_wr=sum(ix_h_wr);
                    
                    phase_fef_cr = []; phase_fef_cr=phase_fef(ix_h_cr,:);
                    phase_it_cr = [];phase_it_cr=phase_it(ix_h_cr,:);
                    diff_phi_cr = [];diff_phi_cr=diff_phi(ix_h_cr,:);
                    rate_it_cr = [];rate_it_cr = rate_it(ix_h_cr,:);
                    val_fef=[];val_it=[]; val_shuff=[];val_shuff=[];
                    [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_cr,phase_it_cr,rate_it_cr,N_sam_sp);
                    
                    %      [coho,scoho,phaso,urcoho,srcoho,ifcoher] = compute_fries_coherence_match(s_it(ix_h_cr,t1_ld<t_h&t_h<t2_ld),rate_it_cr,100,N_sam_sp)
                    
                    val_fef_it_cr(nn,condi,fi)=(val_fef);
                    val_it_it_cr(nn,condi,fi)=(val_it);
                    val_shuff_it_cr(nn,condi,fi)=(val_shuff);
                    
                    if N_wr>0
                        phase_fef_wr = []; phase_fef_wr=phase_fef(ix_h_wr,:);
                        phase_it_wr = [];phase_it_wr=phase_it(ix_h_wr,:);
                        diff_phi_wr = [];diff_phi_wr=diff_phi(ix_h_wr,:);
                        rate_it_wr = [];rate_it_wr = rate_it(ix_h_wr,:);
                        val_fef=[];val_it=[]; val_shuff=[];
                        [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_wr,phase_it_wr,rate_it_wr,N_sam_sp);
                        
                        val_fef_it_wr(nn,condi,fi)=(val_fef);
                        val_it_it_wr(nn,condi,fi)=(val_it);
                        val_shuff_it_wr(nn,condi,fi)=(val_shuff);
                    else
                        val_fef_it_wr(nn,condi,fi)=nan;
                        val_it_it_wr(nn,condi,fi)=nan;
                        val_shuff_it_wr(nn,condi,fi)=nan;
                    end
                    
                    
                end
            end
        end
        Info_wa(ss).val_fef_it_cr=val_fef_it_cr;
        Info_wa(ss).val_it_it_cr=val_it_it_cr;
        Info_wa(ss).val_shuff_it_cr=val_shuff_it_cr;
        
        Info_wa(ss).val_fef_it_wr=val_fef_it_wr;
        Info_wa(ss).val_it_it_wr=val_it_it_wr;
        Info_wa(ss).val_shuff_it_wr=val_shuff_it_wr;
        
        val_shuff_fef_cr = [];
        val_shuff_fef_wr = [];
        
        val_fef_fef_cr = [];
        val_it_fef_cr = [];
        val_diff_fef_cr = [];
        val_fef_fef_n_cr = [];
        val_it_fef_n_cr = [];
        val_diff_fef_n_cr = [];
        
        val_fef_fef_wr = [];
        val_it_fef_wr = [];
        val_diff_fef_wr = [];
        val_fef_fef_n_wr = [];
        val_it_fef_n_wr = [];
        val_diff_fef_n_wr = [];
        obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        FEF=resp(ss).FEF;
        for nn=1:size(FEF,2)
            
            fef_tag=[fef_tag ss*10+nn];
            fef_=FEF{nn}(trii_h,:);
            
            rate_fef= (fef_(:,t1_ld<t_h_sp&t_h_sp<t2_ld));
            n_trials=size(rate_fef,1);
            for fi=1:length(freqs2use)
                
                val_fef_fef_=[];val_it_fef_=[];val_diff_fef_=[];
                phase_fef =squeeze(phase_fef_(:,fi,t1_ld<t_h&t_h<t2_ld));
                phase_it =squeeze(phase_it_(:,fi,t1_ld<t_h&t_h<t2_ld));
                diff_phi=squeeze(diff_phi_(:,fi,t1_ld<t_h&t_h<t2_ld));
                for condi=1:2
                    if condi<2; condi_h=obj_ix_h;else condi_h=obj_ix_h+3; end
                    ix_h_cr = [];
                    ix_h_cr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==1)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    ix_h_wr = [];
                    ix_h_wr = ismember(lfp.condition{:}(trii_h),condi_h)&(lfp.bhv{:}(trii_h)==0)&(~lfp.FEF_artifact(trii_h))&(~lfp.IT_artifact(trii_h));
                    
                    N_cr=sum(ix_h_cr);N_wr=sum(ix_h_wr);
                    
                    phase_fef_cr = []; phase_fef_cr=phase_fef(ix_h_cr,:);
                    phase_it_cr = [];phase_it_cr=phase_it(ix_h_cr,:);
                    diff_phi_cr = [];diff_phi_cr=diff_phi(ix_h_cr,:);
                    rate_fef_cr = [];rate_fef_cr = rate_fef(ix_h_cr,:);
                    val_fef=[];val_it=[];val_diff=[]; val_fef_n=[];val_it_n=[];val_diff_n=[];
                    [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_cr,phase_it_cr,rate_fef_cr,N_sam_sp);
                    
                    val_fef_fef_cr(nn,condi,fi)=(val_fef);
                    val_it_fef_cr(nn,condi,fi)=(val_it);
                    val_shuff_fef_cr(nn,condi,fi)=(val_shuff);
                    
                    if N_wr>0
                        phase_fef_wr = []; phase_fef_wr=phase_fef(ix_h_wr,:);
                        phase_it_wr = [];phase_it_wr=phase_it(ix_h_wr,:);
                        diff_phi_wr = [];diff_phi_wr=diff_phi(ix_h_wr,:);
                        rate_fef_wr = [];rate_fef_wr = rate_fef(ix_h_wr,:);
                        val_fef=[];val_it=[];val_diff=[]; val_fef_n=[];val_it_n=[];val_diff_n=[];
                        [val_fef,val_it,val_shuff]=spl_per_condition_bahmani(phase_fef_wr,phase_it_wr,rate_fef_wr,N_sam_sp);
                        
                        val_fef_fef_wr(nn,condi,fi)=(val_fef);
                        val_it_fef_wr(nn,condi,fi)=(val_it);
                        val_shuff_fef_wr(nn,condi,fi)=(val_shuff);
                        
                    else
                        val_fef_fef_wr(nn,condi,fi)=nan;
                        val_it_fef_wr(nn,condi,fi)=nan;
                        val_shuff_fef_wr(nn,condi,fi)=nan;
                        
                        
                    end
                    
                    
                    
                    
                    
                end
                
            end
        end
        
        Info_wa(ss).val_fef_fef_cr=val_fef_fef_cr;
        Info_wa(ss).val_it_fef_cr=val_it_fef_cr;
        Info_wa(ss).val_shuff_fef_cr=val_shuff_fef_cr;
        
        Info_wa(ss).val_fef_fef_wr=val_fef_fef_wr;
        Info_wa(ss).val_it_fef_wr=val_it_fef_wr;
        Info_wa(ss).val_shuff_fef_wr=val_shuff_fef_wr;
        
        
        toc
    end
    
    
%     save('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\spl_ppc_per_session_bahmani_full_withinarea_pairs.mat'...
%         , 'Info','freqs2use')
end


%% load
% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\si_pval_it.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis\Mat\SPL\spl_ppc_per_session_bahmani_full.mat')
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
for ss=1:size(gen_lfp_pairs,2)
    
    lfp=gen_lfp_pairs(ss);
    %     trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
    
    
    
    IT=[]; IT=resp(ss).IT;
    for nn=1:size(IT,2)
        obj_ix_h=[resp(ss).pref_obj_su(nn) resp(ss).npref_obj_su(nn)];
        obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        
        pp_it=pp_it+1;
        
        spl_it_it_cr(pp_it,:,:)=Info_wa(ss).val_it_it_cr(nn,:,:);
        spl_it_it_wr(pp_it,:,:)=Info_wa(ss).val_it_it_wr(nn,:,:);
        spl_fef_it_cr(pp_it,:,:)=Info_aa(ss).val_fef_it_cr(nn,:,:);
        spl_fef_it_wr(pp_it,:,:)=Info_aa(ss).val_fef_it_wr(nn,:,:);
        %       spl_diff_it_cr(pp_it,:,:)=Info(ss).val_diff_it_cr(nn,:,:);
        %       spl_diff_it_wr(pp_it,:,:)=Info(ss).val_diff_it_wr(nn,:,:);
        %
        spl_it_it_n_cr(pp_it,:,:)=Info_wa(ss).val_it_it_cr(nn,:,:)-Info_wa(ss).val_shuff_it_cr(nn,:,:);
        spl_it_it_n_wr(pp_it,:,:)=Info_wa(ss).val_it_it_wr(nn,:,:)-Info_wa(ss).val_shuff_it_wr(nn,:,:);
        %
        
        spl_fef_it_n_cr(pp_it,:,:)=Info_aa(ss).val_fef_it_cr(nn,:,:)-Info_aa(ss).val_shuff_it_cr(nn,:,:);
        spl_fef_it_n_wr(pp_it,:,:)=Info_aa(ss).val_fef_it_wr(nn,:,:)-Info_aa(ss).val_shuff_it_wr(nn,:,:);
        %       spl_diff_it_n_cr(pp_it,:,:)=Info(ss).val_diff_it_n_cr(nn,:,:);
        %       spl_diff_it_n_wr(pp_it,:,:)=Info(ss).val_diff_it_n_wr(nn,:,:);
        %
        ind_ses_it2(pp_it)=ss;
        
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition,obj_ix_h(1)+[0 3]),:);
        p_val_it(pp_it,1)=signrank(mean(var_h(:,t1:t2),2),mean(var_h(:,1:500),2));
        
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition,obj_ix_h(2)+[0 3]),:);
        p_val_it(pp_it,2)=signrank(mean(var_h(:,t1:t2),2),mean(var_h(:,1:500),2));
        
        var_h=[]; var_h=IT{nn}(ismember(resp(ss).condition,obj_ix_h(3)+[0 3]),:);
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
for ss=1:size(gen_lfp_pairs,2)
    
    lfp=gen_lfp_pairs(ss);
    %     trii_h=[]; trii_h=(~lfp.FEF_artifact)&(~lfp.IT_artifact);
    
    
    obj_ix_h=[lfp.preff_obj lfp.npreff_obj];
    obj_ix_h=[obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
    
    FEF=[]; FEF=resp(ss).FEF;
    for nn=1:size(FEF,2)
        
        pp_fef=pp_fef+1;
        
        spl_it_fef_cr(pp_fef,:,:)=Info_aa(ss).val_it_fef_cr(nn,:,:);
        spl_it_fef_wr(pp_fef,:,:)=Info_aa(ss).val_it_fef_wr(nn,:,:);
        spl_fef_fef_cr(pp_fef,:,:)=Info_wa(ss).val_fef_fef_cr(nn,:,:);
        spl_fef_fef_wr(pp_fef,:,:)=Info_wa(ss).val_fef_fef_wr(nn,:,:);
        %       spl_diff_fef_cr(pp_fef,:,:)=Info(ss).val_diff_fef_cr(nn,:,:);
        %       spl_diff_fef_wr(pp_fef,:,:)=Info(ss).val_diff_fef_wr(nn,:,:);
        
        spl_it_fef_n_cr(pp_fef,:,:)=Info_aa(ss).val_it_fef_cr(nn,:,:)-Info_aa(ss).val_shuff_fef_cr(nn,:,:);
        spl_it_fef_n_wr(pp_fef,:,:)=Info_aa(ss).val_it_fef_wr(nn,:,:)-Info_aa(ss).val_shuff_fef_wr(nn,:,:);
        spl_fef_fef_n_cr(pp_fef,:,:)=Info_wa(ss).val_fef_fef_cr(nn,:,:)-Info_wa(ss).val_shuff_fef_cr(nn,:,:);
        spl_fef_fef_n_wr(pp_fef,:,:)=Info_wa(ss).val_fef_fef_wr(nn,:,:)-Info_wa(ss).val_shuff_fef_wr(nn,:,:);
        %       spl_diff_fef_n_cr(pp_fef,:,:)=Info(ss).val_diff_fef_n_cr(nn,:,:);
        %       spl_diff_fef_n_wr(pp_fef,:,:)=Info(ss).val_diff_fef_n_wr(nn,:,:);
        ind_ses_fef2(pp_fef)=ss;
        
    end
    
end

%% Session selection for comparing corect and Wrong
load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_full.mat');

% load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_500msbaseline_alldata.mat');

psth_cr_it=Val_cr_IT_psth;psth_wr_it=Val_wr_IT_psth;
psth_cr_fef=Val_cr_FEF_psth;psth_wr_fef=Val_wr_FEF_psth;

%%
%status = 1 All , status = 2 SFN; status = 3 performance
ind_ses_itphase=1;%session_selection(4);
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

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef_itphase=setdiff(ixx_fef,ih_s);
fef_itphase_ind_m1=(ind_ses_fef2(ixx_fef_itphase)<52);
fef_itphase_ind_m2=(ind_ses_fef2(ixx_fef_itphase)>51);

%%
ind_ses_fefphase=1;%session_selection(2);   %% FOr preferance put it 7

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
%% Preview niceplots
condi=1;
if noramlized
    %
    bb1_=-0.05;bb2_=0.1;
    
    preview_niceplot_bhv_obj_loc(spl_it_it_n_cr(ixx_it_itphase,:,:),spl_it_it_n_wr(ixx_it_itphase,:,:),freqs2use,'IT unit IT phase ','Norm SPL (a.u)',win_,condi,bb1_,bb2_)
    %
    preview_niceplot_bhv_obj_loc(spl_fef_it_n_cr(ixx_it_fefphase,:,:),spl_fef_it_n_wr(ixx_it_fefphase,:,:),freqs2use,'IT unit FEF phase ','Norm SPL (a.u)',win_,condi,bb1_,bb2_)
    %
    preview_niceplot_bhv_obj_loc(spl_fef_fef_n_cr(ixx_fef_fefphase,:,:),spl_fef_fef_n_wr(ixx_fef_fefphase,:,:),freqs2use,'FEF unit FEF phase ','Norm SPL (a.u)',win_,condi,bb1_,bb2_)
    %
    preview_niceplot_bhv_obj_loc(spl_it_fef_n_cr(ixx_fef_itphase,:,:),spl_it_fef_n_wr(ixx_fef_itphase,:,:),freqs2use,'FEF unit IT phase ','Norm SPL (a.u)',win_,condi,bb1_,bb2_)
    %
    
else
    bb1_=0.21;bb2_=0.31;
    
    
    preview_niceplot_bhv_obj_loc(spl_it_it_cr(ixx_it_itphase,:,:),spl_it_it_wr(ixx_it_itphase,:,:),freqs2use,'IT unit IT phase ',' SPL (a.u)',win_,condi,bb1_,bb2_)
    preview_niceplot_bhv_obj_loc(spl_fef_it_cr(ixx_it_fefphase,:,:),spl_fef_it_wr(ixx_it_fefphase,:,:),freqs2use,'IT unit FEF phase ',' SPL (a.u)',win_,condi,bb1_,bb2_)
    %     % %
    preview_niceplot_bhv_obj_loc(spl_fef_fef_cr(ixx_fef_fefphase,:,:),spl_fef_fef_wr(ixx_fef_fefphase,:,:),freqs2use,'FEF unit FEF phase ',' SPL (a.u)',win_,condi,bb1_,bb2_)
    preview_niceplot_bhv_obj_loc(spl_it_fef_cr(ixx_fef_itphase,:,:),spl_it_fef_wr(ixx_fef_itphase,:,:),freqs2use,'FEF unit IT phase ',' SPL (a.u)',win_,condi,bb1_,bb2_)
end

%% Preview scatters
if noramlized
    bb1_=-0.2;bb2_=0.2;
    %     preview_scatter_bhv_obj_loc(spl_it_fef_n_cr(ixx_fef_itphase,:,:),spl_it_fef_n_wr(ixx_fef_itphase,:,:),freqs2use,'FEF unit IT phase ',bb1_,bb2_,freq_bands,condi)
    
    %     preview_scatter_bhv_obj_loc(spl_it_it_n_cr(ixx_it_itphase,:,:),spl_it_it_n_wr(ixx_it_itphase,:,:),freqs2use,'IT unit IT phase ',bb1_,bb2_,freq_bands,condi)
    
    %     preview_scatter_bhv_obj_loc(spl_fef_it_n_cr(ixx_it_fefphase,:,:),spl_fef_it_n_wr(ixx_it_fefphase,:,:),freqs2use,'IT unit FEF phase ',bb1_,bb2_,freq_bands,condi)
    %
    %     preview_scatter_bhv_obj_loc(spl_fef_fef_n_cr(ixx_fef_fefphase,:,:),spl_fef_fef_n_wr(ixx_fef_fefphase,:,:),freqs2use,'FEF unit FEF phase ',bb1_,bb2_,freq_bands,condi)
    
else
    bb1_=0.1;bb2_=0.4;
        
 %      preview_scatter_bhv_obj_loc(spl_it_it_cr(ixx_it_itphase,:,:),spl_it_it_wr(ixx_it_itphase,:,:),freqs2use,'IT unit IT phase ',bb1_,bb2_,freq_bands,it_itphase_ind_m1,it_itphase_ind_m2)
       preview_scatter_bhv_obj_loc(spl_fef_it_cr(ixx_it_fefphase,:,:),spl_fef_it_wr(ixx_it_fefphase,:,:),freqs2use,'IT unit FEF phase ',bb1_,bb2_,freq_bands,it_fefphase_ind_m1,it_fefphase_ind_m2)
    %
                preview_scatter_bhv_obj_loc(spl_it_fef_cr(ixx_fef_itphase,:,:),spl_it_fef_wr(ixx_fef_itphase,:,:),freqs2use,'FEF unit IT phase ',bb1_,bb2_,freq_bands,fef_itphase_ind_m1,fef_itphase_ind_m2)
       %            preview_scatter_bhv_obj_loc(spl_fef_fef_cr(ixx_fef_fefphase,:,:),spl_fef_fef_wr(ixx_fef_fefphase,:,:),freqs2use,'FEF unit FEF phase ',bb1_,bb2_,freq_bands,fef_fefphase_ind_m1,fef_fefphase_ind_m2)
end

 %% Sample SPL Unit IT phase FEF
% 
% sui=144;
% % for sui=ixx_it
% 
% Color_bhv=[0.7,0.7,0;0,0.7,0.7];
% figure
% hold on
% plot(freqs2use,SmoothN(squeeze(spl_fef_it_cr(sui,1,:)),win_),'color',Color_bhv(1,:),'linewidth',4)
% plot(freqs2use,SmoothN(squeeze(spl_fef_it_wr(sui,1,:)),win_),'color',Color_bhv(2,:),'linewidth',4)
% legend('Correct','Wrong')
% xlabel('Frequency (Hz)')
% ylabel('SPL (a.u)')
% set(gca,'fontsize',14,'fontweight','bold')
% xlim([4 60])
% ylim([0.1 0.35])
% % sui
% % pause
% % close all
% % end
% 
% %% Sample SPL Unit FEF phase IT
% 
% sui=14;
% 
% Color_bhv=[0.7,0.7,0;0,0.7,0.7];
% figure
% hold on
% plot(freqs2use,SmoothN(squeeze(spl_it_fef_cr(sui,1,:)),win_),'color',Color_bhv(1,:),'linewidth',4)
% plot(freqs2use,SmoothN(squeeze(spl_it_fef_wr(sui,1,:)),win_),'color',Color_bhv(2,:),'linewidth',4)
% legend('Correct','Wrong')
% xlabel('Frequency (Hz)')
% ylabel('SPL (a.u)')
% set(gca,'fontsize',14,'fontweight','bold')
% xlim([4 60])
% ylim([0.1 0.35])
% 


end