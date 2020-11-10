function auc= M_index(pref,nonp)

if (pref==0); auc=nan;return;end
if (nonp==0); auc=nan;return;end
auc=(nanmean(pref)-nanmean(nonp))./(nanmean(pref)+nanmean(nonp));

%     auc=0;
% end
end
