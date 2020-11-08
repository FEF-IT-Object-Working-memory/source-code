function resp=LFP_Artifact_remove(resp)
%% Artifact
Th=3;artifacts=[];preview=0;
sigma=5;
gen_lfp_pairs=resp;
gen_lfp_pairs_=resp;
%gen_lfp_pairs.artifact=[];
for ss=1:size(resp,2)
    
    
    st_ind=1;end_ind=2500;
   
        val_it_max=max((gen_lfp_pairs_(ss).LFP_IT(:,st_ind:end_ind)),[],2);
        val_it_min=min((gen_lfp_pairs_(ss).LFP_IT(:,st_ind:end_ind)),[],2);
        val_it=val_it_max-val_it_min;
        th_it=median(val_it)+sigma*std(val_it);
        val_it=max((gen_lfp_pairs_(ss).LFP_IT(:,st_ind:end_ind)),[],2)-min((gen_lfp_pairs_(ss).LFP_IT(:,st_ind:end_ind)),[],2);
        gen_lfp_pairs(ss).IT_artifact=gen_lfp_pairs(ss).bhv;
        gen_lfp_pairs(ss).IT_artifact(abs(val_it)>th_it)=2;
        gen_lfp_pairs(ss).IT_artifact(gen_lfp_pairs(ss).IT_artifact~=2)=0;
        gen_lfp_pairs(ss).IT_artifact(gen_lfp_pairs(ss).IT_artifact==2)=1;
        
   
    
    val_fef_max=max((gen_lfp_pairs_(ss).LFP_FEF(:,st_ind:end_ind)),[],2);
    val_fef_min=min((gen_lfp_pairs_(ss).LFP_FEF(:,st_ind:end_ind)),[],2);
    val_fef=val_fef_max-val_fef_min;
    th_fef=median(val_fef)+sigma*std(val_fef);
    val_fef=max((gen_lfp_pairs_(ss).LFP_FEF(:,st_ind:end_ind)),[],2)-min((gen_lfp_pairs_(ss).LFP_FEF(:,st_ind:end_ind)),[],2);
    
    
    
    gen_lfp_pairs(ss).FEF_artifact=gen_lfp_pairs(ss).bhv;
    gen_lfp_pairs(ss).FEF_artifact(abs(val_fef)> th_fef)=2;
    gen_lfp_pairs(ss).FEF_artifact(gen_lfp_pairs(ss).FEF_artifact~=2)=0;
    gen_lfp_pairs(ss).FEF_artifact(gen_lfp_pairs(ss).FEF_artifact==2)=1;
    
    
    
    if preview
        subplot(613)
        ihx=(gen_lfp_pairs(ss).IT_artifact|gen_lfp_pairs(ss).FEF_artifact);
        plot((gen_lfp_pairs_(ss).FEF(ihx,:))')
        title(strcat('FEF trials=',num2str(size((gen_lfp_pairs_(ss).LFP_FEF(ihx,:)),1))))
        % axis([1 2100 -10 10])
        
        subplot(614)
        if ss<23
            plot((gen_lfp_pairs_(ss).LFP_IT(ihx,:))')
            title(strcat('IT trials=',num2str(size((gen_lfp_pairs_(ss).LFP_IT(ihx,:)),1))))
            %axis([1 2100 -10 10])
            
        end
        ihx=~(gen_lfp_pairs(ss).IT_artifact|gen_lfp_pairs(ss).FEF_artifact);
        subplot(615)
        plot((gen_lfp_pairs_(ss).LFP_FEF(ihx,:))')
        title(strcat('FEF trials=',num2str(size((gen_lfp_pairs_(ss).LFP_FEF(ihx,:)),1))))
        % axis([1 2100 -10 10])
        
        subplot(616)
        plot((gen_lfp_pairs_(ss).IT(ihx,:))')
        title(strcat('IT trials=',num2str(size((gen_lfp_pairs_(ss).LFP_IT(ihx,:)),1))))
        % axis([1 2100 -10 10])
    end
    
    
end

%% remove 50/60 Hz noise

gen_lfp_pairs=noth_filter_lfp(gen_lfp_pairs,size(gen_lfp_pairs,2),size(gen_lfp_pairs,2));

resp=gen_lfp_pairs;
end