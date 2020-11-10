function [test train] = gen_fx_get_equal_part(grp,rate)
catNo = unique(grp);
asiz = []; for cat=catNo'; asiz = [asiz sum(grp == cat)];end
minL = min(asiz);
TrSiz = floor(rate*minL);
TeSiz = minL - TrSiz; 

train = zeros(size(grp)); test = zeros(size(grp));
for cat = catNo'
    ix = grp == cat;
    iix = find(ix);
    Rrp = randperm(length(iix));
    train(iix(Rrp(1:TrSiz))) = 1;
    test(iix(Rrp(end-TeSiz+1:end))) = 1;
end

train = logical(train);
test = logical(test);