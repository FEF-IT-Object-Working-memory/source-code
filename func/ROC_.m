function auc= ROC_(pref,nonp)

if (pref==0); auc=nan;return;end
if (nonp==0); auc=nan;return;end

arr=sort([pref; nonp]);
clear x;
clear y;
for p=1:size(arr)
    sens=size(find(pref>arr(p)),1) ./ ( size(find(pref>arr(p)),1) + size(find(pref<arr(p)),1));
    spec=size(find(nonp<arr(p)),1) ./ ( size(find(nonp<arr(p)),1) + size(find(nonp>arr(p)),1));
    x(p)=1-spec;
    y(p)=sens;
end;
auc=-trapz(x,y);
% if isnan(auc)
%     auc=0;
% end
end
