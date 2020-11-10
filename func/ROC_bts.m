function auc= ROC_bts(pref,nonp)
n_bts=10;
auc_=[];
for iteri=1:n_bts
   pref_=pref(randi(length(pref),length(pref),1));
   npref_=nonp(randi(length(nonp),length(nonp),1));
   auc_(iteri)=ROC_(pref_,npref_);
end
auc=nanmean(auc_);
end