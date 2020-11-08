function Fig2_PPL(resp)
%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig. 2, Fig S3
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


t1_delay= 600; t2_delay= 1200;                                                % Delay period
t1_vis= 50; t2_vis= 350;                                                      % Visual period
t1_fix= -300; t2_fix= 50;                                                     % Fixation period

t = 1:2100;
t = decimate(t,5,'fir');
t_h = t-500;
similar_caxis = 1;bb1 = -1;bb2 = 1;bb1_ = -3; bb2_ = 3;

%% Load
gen_lfp_pairs(1).IT=resp(1).LFP_IT;
gen_lfp_pairs(1).FEF=resp(1).LFP_FEF;
gen_lfp_pairs(1).bhv{1}=resp(1).bhv;gen_lfp_pairs(1).condition{1}=resp(1).condition';
gen_lfp_pairs(1).preff_obj=resp(1).pref_obj_mu;
gen_lfp_pairs(1).npreff_obj=resp(1).npref_obj_mu;
gen_lfp_pairs(1).FEF_artifact=resp(1).FEF_artifact;
gen_lfp_pairs(1).IT_artifact=resp(1).IT_artifact;

%% Core calculation
disp('Wait! Calculating PPL **********')
%% Core calculation Fig S3 a [Trial matched]
% method trail matching minus 1000 time  shuffle
cal=1;
if cal
    N_shuff=11;
    st_ind=1; end_ind= 2100;
    ind_b=[40 90];
    %freqs2use  =logspace(log10(4),log10(100),60); % 4-100 Hz in 30 steps
    
    % freqs2use  =[1:35 37:2:100]; % 4-100 Hz in 30 steps
    freqs2use  =[4:1:70 71:5:130];                                         % frequency range
    
    %   load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
    %
    ppl_cr=[];ppl_wr=[];plv_cr_z=[];plv_wr_z=[];
    ppl_cr_n=[];ppl_wr_n=[];plv_cr_n_z=[];plv_wr_n_z=[];
    plv_cr_shuff_m=[];plv_wr_shuff_m=[];plv_cr_shuff_m_z=[];plv_wr_shuff_m_z=[];
    
    for ss=1:size(gen_lfp_pairs,2)
        
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
        
        
        %% PLV
        %%% Correct
        plv_cr_h=[];
        plv_cr_z_h =[];
        plv_cr_n_h=[];
        plv_cr_n_z_h=[];
        plv_cr_shuff_m_h=[];
        plv_cr_shuff_m_z_h=[];
        
        plv_wr_h=[];
        plv_wr_z_h =[];
        plv_wr_n_h=[];
        plv_wr_n_z_h=[];
        plv_wr_shuff_m_h=[];
        plv_wr_shuff_m_z_h=[];
        
        s_it_h = []; s_it_h = lfp.IT(:,st_ind:end_ind);
        s_fef_h = []; s_fef_h = lfp.FEF(:,st_ind:end_ind);
        P_it = []; P_fef = [];
        [P_it,P_fef,freqs2use_,t]= wavelet_core_DS_pertrial(s_it_h,s_fef_h,freqs2use);
        
        for condi=1:6
            
            N_trials_cr = [];N_trials_cr =  sum(logical(ix_h_cr(condi,:)));
            N_trials_wr = [];N_trials_wr =  sum(logical(ix_h_wr(condi,:)));
            
            if N_trials_cr>3
                phi_it_cr = [];phi_it_cr = squeeze(angle(P_it(logical(ix_h_cr(condi,:)),:,:)));
                phi_fef_cr = [];phi_fef_cr = squeeze(angle(P_fef(logical(ix_h_cr(condi,:)),:,:)));
            end
            if N_trials_wr>3
                phi_it_wr = [];phi_it_wr = squeeze(angle(P_it(logical(ix_h_wr(condi,:)),:,:)));
                phi_fef_wr = [];phi_fef_wr = squeeze(angle(P_fef(logical(ix_h_wr(condi,:)),:,:)));
                
            end
            
            N_sam=min(N_trials_cr,N_trials_wr);
            N_sam(N_sam<4)=4;
            plv_h_cr_=[];  plv_h_shuff_cr=[];
            plv_h_wr_=[];  plv_h_shuff_wr=[];
            for iteri = 1:N_shuff
                if N_trials_cr>3
                    ind_h = [];ind_h = randperm(N_trials_cr,N_sam);
                    plv_h_cr_(iteri,:,:)=squeeze(abs(nanmean(exp(1j*((phi_fef_cr(ind_h,:,:)-phi_it_cr(ind_h,:,:)))),1)));
                    
                    plv_h_shuff_cr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_cr(ind_h(randperm(N_sam)),:,:)...
                        -phi_it_cr(ind_h(randperm(N_sam)),:,:)))),1)));
                end
                if N_trials_wr>3
                    ind_h = [];ind_h = randperm(N_trials_wr,N_sam);
                    plv_h_wr_(iteri,:,:)=squeeze(abs(nanmean(exp(1j*((phi_fef_wr(ind_h,:,:)-phi_it_wr(ind_h,:,:)))),1)));
                    
                    plv_h_shuff_wr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_wr(ind_h(randperm(N_sam)),:,:)...
                        -phi_it_wr(ind_h(randperm(N_sam)),:,:)))),1)));
                end
                
            end
            
            %%% Fill variables
            if N_trials_cr>3
                plv_cr_h(condi,:,:)=squeeze(nanmean(plv_h_cr_));
                val_ = []; val_ =squeeze(nanmean(plv_h_cr_(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(plv_h_cr_,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(plv_h_cr_,3));
                plv_cr_z_h(condi,:,:)=(squeeze(nanmean(plv_h_cr_))-val_m)./val_std;
                plv_cr_n_h(condi,:,:)=(squeeze(nanmean(plv_h_cr_-plv_h_shuff_cr)));
                var_h=[];var_h=(squeeze((plv_h_cr_-plv_h_shuff_cr)));
                val_ = []; val_ =squeeze(nanmean(var_h(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(var_h,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(var_h,3));
                plv_cr_n_z_h(condi,:,:)=(squeeze(nanmean(var_h))-val_m)./val_std;
                plv_cr_shuff_m_h(condi,:,:)=squeeze(nanmean(plv_h_shuff_cr));
                val_ = []; val_ =squeeze(nanmean(plv_h_shuff_cr(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(plv_h_shuff_cr,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(plv_h_shuff_cr,3));
                plv_cr_shuff_m_z_h(condi,:,:)=(squeeze(nanmean(plv_h_shuff_cr))-val_m)./val_std;
            else
                plv_cr_h(condi,:,:)=nan*ones(67,420);
                plv_cr_z_h(condi,:,:)=nan*ones(67,420);
                plv_cr_n_h(condi,:,:)=nan*ones(67,420);
                plv_cr_n_z_h(condi,:,:)=nan*ones(67,420);
                plv_cr_shuff_m_h(condi,:,:)=nan*ones(67,420);
                plv_cr_shuff_m_z_h(condi,:,:)=nan*ones(67,420);
                
            end
            
            if N_trials_wr>3
                plv_wr_h(condi,:,:)=squeeze(nanmean(plv_h_wr_));
                val_ = []; val_ =squeeze(nanmean(plv_h_wr_(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(plv_h_wr_,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(plv_h_wr_,3));
                plv_wr_z_h(condi,:,:)=(squeeze(nanmean(plv_h_wr_))-val_m)./val_std;
                plv_wr_n_h(condi,:,:)=(squeeze(nanmean(plv_h_wr_-plv_h_shuff_wr)));
                var_h=[]; var_h=(squeeze((plv_h_wr_-plv_h_shuff_wr)));
                val_ = []; val_ =squeeze(nanmean(var_h(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(var_h,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(var_h,3));
                plv_wr_n_z_h(condi,:,:)=(squeeze(nanmean(var_h))-val_m)./val_std;
                plv_wr_shuff_m_h(condi,:,:)=squeeze(nanmean(plv_h_shuff_wr));
                val_ = []; val_ =squeeze(nanmean(plv_h_shuff_wr(:,:,ind_b(1):ind_b(2)),3));
                val_m=nanmean(val_,1);val_m=repmat(val_m',1,size(plv_h_shuff_wr,3));
                val_std=nanstd(val_,[],1);val_std=repmat(val_std',1,size(plv_h_shuff_wr,3));
                plv_wr_shuff_m_z_h(condi,:,:)=(squeeze(nanmean(plv_h_shuff_wr))-val_m)./val_std;
                
            else
                plv_wr_h(condi,:,:)=nan*ones(67,420);
                plv_wr_z_h(condi,:,:)=nan*ones(67,420);
                plv_wr_n_h(condi,:,:)=nan*ones(67,420);
                plv_wr_n_z_h(condi,:,:)=nan*ones(67,420);
                plv_wr_shuff_m_h(condi,:,:)=nan*ones(67,420);
                plv_wr_shuff_m_z_h(condi,:,:)=nan*ones(67,420);
                
            end
            
        end
        
        ppl_cr(ss,:,:,:)=plv_cr_h;
        plv_cr_z(ss,:,:,:)=plv_cr_z_h;
        ppl_cr_n(ss,:,:,:)=plv_cr_n_h;
        plv_cr_n_z(ss,:,:,:)=plv_cr_n_z_h;
        plv_cr_shuff_m(ss,:,:,:)=plv_cr_shuff_m_h;
        plv_cr_shuff_m_z(ss,:,:,:)=plv_cr_shuff_m_z_h;
        
        ppl_wr(ss,:,:,:)=plv_wr_h;
        plv_wr_z(ss,:,:,:)=plv_wr_z_h;
        ppl_wr_n(ss,:,:,:)=plv_wr_n_h;
        plv_wr_n_z(ss,:,:,:)=plv_wr_n_z_h;
        plv_wr_shuff_m(ss,:,:,:)=plv_wr_shuff_m_h;
        plv_wr_shuff_m_z(ss,:,:,:)=plv_wr_shuff_m_z_h;
        
        
        toc
        
    end
    
    %     save('plv_wavelet_over_trials_shuffle_corrected_trail_matched.mat','t','freqs2use',...
    %         'plv_cr','plv_wr','plv_cr_z','plv_wr_z',...
    %         'plv_cr_n','plv_wr_n','plv_cr_n_z','plv_wr_n_z',...
    %         'plv_cr_shuff_m','plv_wr_shuff_m','plv_cr_shuff_m_z','plv_wr_shuff_m_z')
    % %
    %      save('plv_wavelet_over_trials_shuffle_corrected_trail_matched_full.mat','t','freqs2use',...
    %         'plv_cr','plv_wr',...
    %         'plv_cr_n','plv_wr_n',...
    %         'plv_cr_shuff_m','plv_wr_shuff_m')
    
end
ppl_cr_tri_matched = ppl_cr_n;ppl_wr_tri_matched = ppl_wr_n;

%%% shuffle within conditions
cal = 1;
if cal
    %     ind_b = [1 110];                                                         % base line when down sample fs to 200 ( -300 -50)ms
    %     load gen_lfp_pairs_500ms_baseline_final                                       % load lfp data
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')
    %     load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_lfp_pairs_500ms_baseline_full.mat')
    
    
    
    N_shuff = 11;
    
    freqs2use   = [4:1:70 71:5:130];                                         % frequency range
    st_ind = 1; end_ind =  2100;                                               % time of task when fixatio is 500 ms
    
    %%% Definition of variabels
    ppl_cr = [];ppl_wr = [];                                                   % raw PLV for Cr& Wr
    plv_cr_z = [];plv_wr_z = [];                                               % baseline nomalized raw PLV for Cr& Wr
    ppl_cr_n = [];ppl_wr_n = [];                                               % shuffle corrected PLV for Cr& Wr
    plv_cr_n_z = [];plv_wr_n_z = [];                                           % shuffle corrected baseline nomalized PLV for Cr& Wr
    plv_shuff_cr_m = [];plv_shuff_wr_m = [];                                   % mean of shuffled  PLV for Cr& Wr
    plv_shuff_cr_std = [];plv_shuff_wr_std = [];                               % mean of shuffled baesline normalied PLV for Cr& Wr
    NS = size(gen_lfp_pairs,2);
    for ss = 1:size(gen_lfp_pairs,2)
        
        
        tic
        %%% Determining  Condtions index
        lfp = gen_lfp_pairs(ss);
        ix_h_cr = []; ix_h_wr = [];
        obj_ix_h = [lfp.preff_obj lfp.npreff_obj];
        obj_ix_h = [obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2) ];
        for condi=1:6
            if condi<4; condi_h=obj_ix_h(condi);else condi_h=obj_ix_h(condi-3)+3; end
            ix_h_cr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==1)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            ix_h_wr(condi,:)=ismember(lfp.condition{:},condi_h)&(lfp.bhv{:}==0)&(~lfp.FEF_artifact)&(~lfp.IT_artifact);
            
        end
        
        
        %%% calculate morlet wavelt transform and spaontanous phase values
        %%% for each channel
        
        s_it_h = []; s_it_h = lfp.IT(:,st_ind:end_ind);
        s_fef_h = []; s_fef_h = lfp.FEF(:,st_ind:end_ind);
        p_it = []; p_fef = [];
        [p_it,p_fef,freqs2use_,t]= wavelet_core_DS_pertrial(s_it_h,s_fef_h,freqs2use);
        phi_it = [];phi_it=angle(p_it);
        phi_fef =[];phi_fef=angle(p_fef);
        N_trials=size(phi_it,1);
        
        %%% calculate PLV
        
        plv_cr_h=[];
        plv_cr_n_h=[];
        plv_shuff_cr_m_h=[];
        plv_shuff_cr_std_h=[];
        plv_wr_h=[];
        plv_wr_n_h=[];
        plv_shuff_wr_m_h=[];
        plv_shuff_wr_std_h=[];
        
        for condi=1:6
            N_trials_cr = [];N_trials_cr =  sum(logical(ix_h_cr(condi,:)));
            N_trials_wr = [];N_trials_wr =  sum(logical(ix_h_wr(condi,:)));
            
            %             N_sam=min(N_trials_cr,N_trials_wr);
            %             if(N_sam<3);N_sam=3;end
            N_sam=3; %no trial matching
            phi_it_cr = [];phi_it_cr = squeeze(phi_it(logical(ix_h_cr(condi,:)),:,:));
            phi_fef_cr = [];phi_fef_cr = squeeze(phi_fef(logical(ix_h_cr(condi,:)),:,:));
            
            if N_trials_wr>0
                phi_it_wr = [];phi_it_wr = squeeze(phi_it(logical(ix_h_wr(condi,:)),:,:));
                phi_fef_wr = [];phi_fef_wr = squeeze(phi_fef(logical(ix_h_wr(condi,:)),:,:));
            end
            
            %%% Core PLV calculation over trials
            plv_h_cr=[]; plv_h_wr=[];
            
            plv_h_cr=squeeze(abs(nanmean(exp(1j*((phi_fef_cr-phi_it_cr))),1)));
            
            if N_trials_wr>=N_sam
                plv_h_wr=squeeze(abs(nanmean(exp(1j*((phi_fef_wr-phi_it_wr))),1)));
                
            else
                plv_h_wr=nan*ones(size(freqs2use,2),length(t));
            end
            
            %             plv_h_shuff_cr=[];plv_h_shuff_wr=[];
            %             for iteri = 1:N_shuff
            %
            %                 plv_h_shuff_cr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_cr(randperm(N_trials_cr,N_trials_cr),:,:)...
            %                     -phi_it_cr(randperm(N_trials_cr,N_trials_cr),:,:)))),1)));
            %
            %                 if N_trials_wr>=N_sam
            %                     plv_h_shuff_wr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_wr(randperm(N_trials_wr,N_trials_wr),:,:)...
            %                         -phi_it_wr(randperm(N_trials_wr,N_trials_wr),:,:)))),1)));
            %                 else
            %                     plv_h_shuff_wr(iteri,:,:) = nan*ones(size(freqs2use,2),length(t));
            %                 end
            %             end
            %
            plv_h_shuff_cr=[];plv_h_shuff_wr=[];
            if N_trials_wr>=N_sam
                phi_it_shuff=[phi_it_cr;phi_it_wr];phi_fef_shuff=[phi_fef_cr;phi_fef_wr];
            else
                phi_it_shuff=phi_it_cr;phi_fef_shuff=phi_fef_cr;
                
            end
            N_trials_shuff=size(phi_it_shuff,1);
            for iteri = 1:N_shuff
                
                plv_h_shuff_cr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_shuff(randperm(N_trials_shuff,N_trials_cr),:,:)...
                    -phi_it_shuff(randperm(N_trials_shuff,N_trials_cr),:,:)))),1)));
                
                if N_trials_wr>=N_sam
                    plv_h_shuff_wr(iteri,:,:) = squeeze(abs(nanmean(exp(1j*((phi_fef_shuff(randperm(N_trials_shuff,N_trials_wr),:,:)...
                        -phi_it_shuff(randperm(N_trials_shuff,N_trials_wr),:,:)))),1)));
                else
                    plv_h_shuff_wr(iteri,:,:) = nan*ones(size(freqs2use,2),length(t));
                end
            end
            
            %%% fill temporary variables and baseline normalization
            
            plv_cr_h(condi,:,:) = plv_h_cr;
            plv_cr_n_h(condi,:,:)=(plv_h_cr-squeeze(nanmean(plv_h_shuff_cr)));
            plv_shuff_cr_m_h(condi,:,:)=(squeeze(nanmean(plv_h_shuff_cr)));
            plv_shuff_cr_std_h(condi,:,:)=(squeeze(nanstd(plv_h_shuff_cr)));
            
            
            plv_wr_h(condi,:,:) = plv_h_wr;
            plv_wr_n_h(condi,:,:)=(plv_h_wr-squeeze(nanmean(plv_h_shuff_wr)));
            plv_shuff_wr_m_h(condi,:,:)=(squeeze(nanmean(plv_h_shuff_wr)));
            plv_shuff_wr_std_h(condi,:,:)=(squeeze(nanstd(plv_h_shuff_wr)));
            
        end
        
        %%% fill variables
        ppl_cr(ss,:,:,:)=plv_cr_h;
        ppl_cr_n(ss,:,:,:)=plv_cr_n_h;
        plv_shuff_cr_m(ss,:,:,:)=plv_shuff_cr_m_h;
        plv_shuff_cr_std(ss,:,:,:)=plv_shuff_cr_std_h;
        ppl_wr(ss,:,:,:)=plv_wr_h;
        ppl_wr_n(ss,:,:,:)=plv_wr_n_h;
        plv_shuff_wr_m(ss,:,:,:)=plv_shuff_wr_m_h;
        plv_shuff_wr_std(ss,:,:,:)=plv_shuff_wr_std_h;
        toc
        
    end
    
    %     save(strcat('G:\Ehsan\FEF-IT\Mat\PLV\plv_wavelet_over_trials_shuffle_corrected_withincond_final_end_both2020.mat'),'freqs2use',...
    %         'plv_cr','plv_wr','plv_cr_n','plv_wr_n',...
    %         'plv_shuff_cr_m','plv_shuff_cr_std','plv_shuff_wr_m','plv_shuff_wr_std')
    %
    %     save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\PPL_full_shuff_1001.mat'),'freqs2use',...
    %         'plv_cr','plv_wr','plv_cr_n','plv_wr_n',...
    %         'plv_shuff_cr_m','plv_shuff_cr_std','plv_shuff_wr_m','plv_shuff_wr_std')
    %
    %     save(strcat('F:\Data\Reaserch\Thesis\Data Analysis\Mat\PLV\plv_wavelet_over_2020_full_bothdata_11.mat'),'freqs2use',...
    %         'plv_cr','plv_wr','plv_cr_n','plv_wr_n',...
    %         'plv_shuff_cr_m','plv_shuff_cr_std','plv_shuff_wr_m','plv_shuff_wr_std')
end

%% Load data
% load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\PPL\PPL_shuff_1001_trial_matched.mat')

% load(strcat('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Mat\PPL\PPL_shuff_1001_2.mat'));

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
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_wr = (ppl_wr_h-val_b)./val_sd;
end

% Baseline Correction for not shuffle corretion data Fig. S3a
if baseline_corrected
    ppl_cr_h = ppl_cr_nsc;ppl_wr_h = ppl_wr_nsc;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr_nsc = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_wr_nsc = (ppl_wr_h-val_b)./val_sd;
end

% Baseline Correction for Trial matched data Fig. S3a
if baseline_corrected
    ppl_cr_h = ppl_cr_tri_matched;ppl_wr_h = ppl_wr_tri_matched;
    val_b = nanmean(ppl_cr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
    val_b = repmat(val_b,1,1,1,420);
    ppl_cr_tri_matched = (ppl_cr_h-val_b)./val_sd;
    
    val_b = nanmean(ppl_wr_h(:,:,:,ind_b(1):ind_b(2)),4);
    val_sd = nanstd(val_b);
    val_sd = repmat(val_sd,NS,1,1);
    val_sd = repmat(val_sd,1,6,1,420);
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
ind_ses=1;%session_selection(4);                                                               % All sessions  86
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. 2b Beta band PPL during sample and delay differed between Cr and Wr trials
try
ind_ses=1;%session_selection(4);                                                               % All sessions  86
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
catch
end
%% Fig. 2c Comparison of beta band PPL during the delay period for correct versus wrong trials
ind_ses=1;%session_selection(4);                                                               % All sessions  86
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
ind_ses=1;%session_selection(4);                                                               % All sessions  86
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
ind_ses=1;%session_selection(4);                                                               % All sessions  86
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

%% Fig. S3a Beta band PPL was greater for correct trials without the shuffle correction
ind_ses=1;%session_selection(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3a, Rezayat.E, et al')
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
%% Fig. S3a Beta band PPL was greater for correct trials without the shuffle correction
ind_ses=1;%session_selection(4);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3a, Rezayat.E, et al')
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

%% Fig. S3b The trend toward higher beta band PPL on correct trials was present in both M1 and M2
figure('name','Fig. S3b, Rezayat.E, et al')
for iteri=1:2
    ax= subplot(1,2,iteri);
    str=[];
    ix_h = []; if ismember(iteri,[1 3]);ix_h = ind_ses_m1;str='M1'; else ix_h = ind_ses_m2; str='M2'; end
    P = []; P = ppl_cr_smoothed-ppl_wr_smoothed;str=[str ''];
    if size(P,1)>1
        var_h=squeeze((P(ind_ses,1,:,:)));
        var_h=squeeze(nanmean((var_h),1));
    else
        var_h=squeeze((P(ind_ses,1,:,:)));
        var_h=squeeze(((var_h)));
    end
    h=pcolor(t_h/1000,freqs2use,var_h);
    
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

%% Fig. S3c Inter-areal beta PPL encoded object identity during the delay period
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3c, Rezayat.E, et al')
for iteri=1:3
    ax= subplot(1,3,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (ppl_cr_smoothed(:,1,:,:));str = ' Pref';
        case 2
            P = (ppl_cr_smoothed(:,3,:,:));str = ' NPref';
        case 3
            P =( ppl_cr_smoothed(:,1,:,:))-(ppl_cr_smoothed(:,3,:,:));str = ' Pref - NPref';
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S3c Inter-areal beta PPL encoded object identity during the delay period
try
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3c, Rezayat.E, et al')
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
catch
end
%% Fig. S3c Inter-areal beta PPL encoded object identity during the delay period
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
f1=freq_bands(1,3);f2=freq_bands(2,3);
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,3,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. S3c, Rezayat.E, et al')
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

%% Fig. S3d Inter-areal beta PPL encoded the location of the sample during the delay period
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3d, Rezayat.E, et al')
for iteri=1:3
    ax= subplot(1,3,iteri);
    str=[];P = [];
    switch iteri
        case 1
            P = (ppl_cr_smoothed(:,1,:,:));str = ' In';
        case 2
            P = (ppl_cr_smoothed(:,4,:,:));str = ' Out';
        case 3
            P =( ppl_cr_smoothed(:,1,:,:))-(ppl_cr_smoothed(ind_ses,4,:,:));str = ' In - Out';
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
    if iteri==1
        ylabel('Frequency (Hz)')
        xlabel('Time from sample onset (sec.)')
        ax.YTick=[4 15 35 100];
        text(0.5,5,strcat('n=',num2str(length(ind_ses))),'fontsize',14)
    end
    set(gca,'fontsize',14,'fontweight','bold')
end

%% Fig. S3d Inter-areal beta PPL encoded the location of the sample during the delay period
try
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S3c, Rezayat.E, et al')
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
catch
end
%% Fig. S3d Inter-areal beta PPL encoded the location of the sample during the delay period
ind_ses=1;%session_selection(3);                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
f1=freq_bands(1,3);f2=freq_bands(2,3);
val_cr=[];val_cr=squeeze(nanmean(nanmean(ppl_cr(:,1,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
val_wr=[];val_wr=squeeze(nanmean(nanmean(ppl_cr(:,4,freqs2use>f1&freqs2use<f2,t1_delay<t_h&t_h<t2_delay),4),3));
figure('name','Fig. S3c, Rezayat.E, et al')
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

%% State for beta-band PPL during the delay period compared to baseline For Correct trials
try
    disp('%%% State for beta-band PPL during the delay period compared to baseline For Correct trials %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end
%% State for beta-band PPL during the delay period compared to baseline For Wrong trials
try
    disp('%%% State for beta-band PPL during the delay period compared to baseline For Wrong trials %%%')
    ind_ses=session_selection(3);                                               % All sessions  86
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
catch
end
%% State for beta-band PPL during the delay period for Correct trials
try
    disp('%%% State for beta-band PPL during the delay period for Correct trials %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end
%% State for object selectivity of beta-band PPL during the delay period
try
    disp('%%% State for object selectivity of beta-band PPL during the Delay period %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end
%% State for location selectivity of beta-band PPL during the delay period
try
    disp('%%% State for location selectivity of beta-band PPL during the Delay period %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end
%% State for Behavioral correlation of beta-band PPL during the delay period
try
    disp('%%% State for Behavioral correlation of beta-band PPL during the Delay period %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end
%% State for Behavioral correlation of object selectivity during the Delay period
try
    disp('%%% State for Behavioral correlation of object selectivity during the Delay period %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end
%% State for Behavioral correlation of location selectivity during the Delay period
try
    disp('%%% State for Behavioral correlation of location selectivity during the Delay period %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end
%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Trial matched]
try
    disp('%%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Trial matched] %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end
%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Not Shuffle Corrected]
try
    disp('%%% State for Control of Behavioral correlation of beta-band PPL during the delay period [Not Shuffle Corrected] %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end

%% State for beta-band PPL during the visual period for Correct trials
try
    disp('%%% State for beta-band PPL during the visual period for Correct trials %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end

%% State for object selectivity of beta-band PPL during the visual period
try
    disp('%%% State for object selectivity of beta-band PPL during the visual period %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end

%% State for location selectivity of beta-band PPL during the visual period
try
    disp('%%% State for location selectivity of beta-band PPL during the visual period %%%')
    ind_ses=session_selection(3);                                                               % All sessions  86
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
catch
end

%% State for Behavioral correlation of beta-band PPL during the visual period
try
    disp('%%% State for Behavioral correlation of beta-band PPL during the visual period %%%')
    ind_ses=session_selection(4);                                                               % All sessions  86
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
catch
end