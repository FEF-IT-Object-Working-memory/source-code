function obj_ind=find_preferance(resp_unit,condition)
%% unit object preferance based on Dicarlo paper 
% Fast Readout of Object Identity from Macaque Inferior Temporal Cortex
% Science  04 Nov 2005
win_ = 25; psth =@(x)SmoothN(1000*mean(x,1), win_);
Condition = [[1,2,3];[4,5,6]];
 Neuron_var=[];test_=[];
  var_h = [];IT_neurons=resp_unit;
    for c=1:3; 
        ix_h = []; ix_h = ismember(condition,[c c+3]);
        var_h(c,:)=psth( IT_neurons(ix_h,300:1800));
    end
    var_h_ = var_h-repmat(mean2(var_h(:,1:300)),size(var_h,1),size(var_h,2));
    
     fir_drive_h = []; fir_drive_h = mean(var_h_);
    
    [peak_h peak_ix] = max(fir_drive_h(300:1000));
    peak_ix = 300+peak_ix;
    treh_h = 0.2*peak_h;
    var_ix_h = (fir_drive_h>treh_h);
    var_ix_h_idx_ = [var_ix_h 0; 0 var_ix_h];
    var_ix_h_idx = 2*var_ix_h_idx_(1,:)-var_ix_h_idx_(2,:);
    idx_st = find(var_ix_h_idx==2);idx_end=find(var_ix_h_idx==-1);
    idx_st_h = idx_st(peak_ix-idx_st>-1);idx_h =[]; [~,idx_h]=min(peak_ix-idx_st_h);
    win_st = idx_st_h(idx_h); win_ix_h=find( win_st==idx_st);
    win_end = idx_end(win_ix_h);
    if length(idx_st)>win_ix_h
        while ((win_end-idx_st(win_ix_h+1))<25) & ((win_end-win_st)<475)
            if (idx_end(win_ix_h+1)-win_st)<375
                win_end=idx_end(win_ix_h+1);
                win_ix_h=win_ix_h+1;
                
                if win_ix_h+1> length(idx_st)
                    break;
                end
            else
                break;
            end
        end
    end
    win_st = win_st+300;
    win_end = win_end+300;
Neuron_var=[];obj_ind=[];
    for c=1:3;
        ix_h = []; ix_h = ismember(condition,[c c+3]);
        Neuron_var(c)=mean2( IT_neurons(ix_h,win_st:win_end));
    end
     [~,obj_ind.pref]=max(Neuron_var);
     [~,obj_ind.npref]=min(Neuron_var);  
   
    
end