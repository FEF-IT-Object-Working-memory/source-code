function out = gen_fx_get_crosscorr(I1,I2)
C_cross=xcorr(I1,I2,'coeff');%C_cross=C_cross(20000:end-20000);
C1=xcorr(I1,'coeff');%C1=C1(20000:end-20000);
C2=xcorr(I2,'coeff');%C2=C2(20000:end-20000);
out=max(C_cross);


% Pt = [];
% Tu = [];
% thr = 10000;
% for cnt = 1:rep
%     if length(grp0) > thr
%         ix = randperm(length(grp0),thr);
%         grp = grp0(ix); I = I0(ix,:);
%     else
%         grp = grp0; I = I0;
%     end
%     
%     [test train] = gen_fx_get_equal_part(grp,rate);
%     if sum(test)>5 && sum(train) > 5
%         cls = gen_fx_MC_SVM(I(test,:),I(train,:),grp(train));
%         Pt = [Pt; sum(cls == grp(test))/sum(test)];
%         Ct(:,:,cnt)= confusionmat(cls,grp(test));
%         Tu = [Tu diag(Ct(:,:,cnt)) ./ (sum(Ct(:,:,cnt)))'];
%     else
%         Pt = [Pt; nan];
%         Ct(:,:,cnt)= nan;
%         Tu = [Tu nan];
%     end
% end
% out.C = Ct; out.pt = Pt; out.tu = Tu;
end