% sorting criterial

%%
clear all
rep = 30;
rate = .7;
win_siz = 300; % size of spikes in one window

bin_no = 50;
%%
kopol_path='F:\Data\Reaserch\Thesis\FEF_IT project\Kopol data\Spike Waveform\';
folder_names= dir(kopol_path);
folder_names=folder_names(3:end);

load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')
load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_resp_pairs_500ms_baseline_full.mat')

%%
%
load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\Normalized_firing_rate_full.mat');

load ('F:\Data\Reaserch\Thesis\Data Analysis\Mat\Firing Rate\firing_rate_full.mat');

psth_cr_it_n=Val_cr_IT_psth_norm;psth_wr_it_n=Val_wr_IT_psth_norm;
psth_cr_fef_n=Val_cr_FEF_psth_norm;psth_wr_fef_n=Val_wr_FEF_psth_norm;
%
% psth_cr_it_n=Val_cr_IT_psth;psth_wr_it_n=Val_wr_IT_psth;
% psth_cr_fef_n=Val_cr_FEF_psth;psth_wr_fef_n=Val_wr_FEF_psth;

psth_cr_it=Val_cr_IT_psth;psth_wr_it=Val_wr_IT_psth;
psth_cr_fef=Val_cr_FEF_psth;psth_wr_fef=Val_wr_FEF_psth;


th_h=0.002;p_val_d=0.05;
ixx_it=1:253;ixx_fef=1:180;

ih_s=ixx_it(find(abs(mean(squeeze(psth_wr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);
ih_s=ixx_it(find(abs(mean(squeeze(psth_cr_it(ixx_it,1,500:1800)),2))<=th_h));
ixx_it=setdiff(ixx_it,ih_s);

ih_s=ixx_fef(find(abs(mean(squeeze(psth_wr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ih_s=ixx_fef(find(abs(mean(squeeze(psth_cr_fef(ixx_fef,1,500:1800)),2))<=th_h));
ixx_fef=setdiff(ixx_fef,ih_s);
ind_ses_it_included=zeros(253,1);ind_ses_it_included(ixx_it)=1;
ind_ses_fef_included=zeros(180,1);ind_ses_fef_included(ixx_fef)=1;
it_uni_inc=[];fef_uni_inc=[];
for ss=1:92
   it_uni_inc{ss}=ones(sum(ismember(ind_ses_it,ss)),1);
   fef_uni_inc{ss}=ones(sum(ismember(ind_ses_fef,ss)),1);

end
%%
cnt_it=0;cnt_fef=0;
for ncnt=1:86%size(gen_session_inf,2) % this loop works for all sessions
   
    ncnt
    tic1=tic;
    if ncnt<52
    path_h=strcat(kopol_path,gen_session_inf(ncnt).name,'\');
    else
       stanford_path= 'F:\Data\Reaserch\Thesis\FEF_IT project\Sanford data\Spike waveforms';
        path_h=strcat(stanford_path,gen_session_inf(ncnt).name,'\');
    end
    
    %% IT
    load(strcat(path_h,'IT_clsInf_plx.mat'))
    load(strcat(path_h,'SP_IT_plx.mat'))
    grp = clsInf.cls;
    cls = clsInf.Units;
%     cnt_it=cnt_it+length(cls);
    for cc = 1:length(cls); IT_len{ncnt}(cc) = sum(grp==cls(cc)); end
    if length(cls) > 2
        cls(cls==0) = [];
        ix = (grp ~= 0);
    else
        ix = (grp > -1);
    end
    
    rc = nchoosek(1:length(cls),2);
    for pcnt=1:size(rc,1)
        ix0=[]; ix0 = ix & ismember(grp,cls(rc(pcnt,:)));
        I = []; I = SP(ix0,:); g1=[]; g1 = grp(ix0);
        t_itT(ncnt,pcnt) = gen_fx_get_svm_sor(g1,I,rate,rep);
        
        tim_bin = round(linspace(1,size(I,1)-1,bin_no));
        for tcnt=1:(length(tim_bin)-1) 
            ix2=[]; ix2 = tim_bin(tcnt):tim_bin(tcnt+1);
            t_it(ncnt,pcnt,tcnt) = gen_fx_get_svm_sor(g1(ix2),I(ix2,:),rate,rep);
        end
    end
    
    
    %% FEF
     load(strcat(path_h,'FEF_clsInf_plx.mat'))
    load(strcat(path_h,'SP_FEF_plx.mat'))
    grp = clsInf.cls;
    cls = clsInf.Units;
    for cc = 1:length(cls); FEF_len{ncnt}(cc) = sum(grp==cls(cc)); end
%     if length(cls) > 2
%         cls(cls==0) = [];
%         ix = (grp ~= 0);
%     else
        ix = (grp > -1);
%     end
    rc = nchoosek(1:length(cls),2);
    for pcnt=1:size(rc,1)
        ix0=[]; ix0 = ix & ismember(grp,cls(rc(pcnt,:)));
        I = []; I = SP(ix0,:); g1=[]; g1 = grp(ix0);
        t_fefT(ncnt,pcnt) = gen_fx_get_svm_sor(g1,I,rate,rep);
        
        tim_bin = round(linspace(1,size(I,1)-1,bin_no));
        for tcnt=1:(length(tim_bin)-1) 
            ix2=[]; ix2 = tim_bin(tcnt):tim_bin(tcnt+1);
            t_fef(ncnt,pcnt,tcnt) = gen_fx_get_svm_sor(g1(ix2),I(ix2,:),rate,rep);
        end
    end
    
    toc(tic1)
end