% size (trial X time)


function niceplot(input,t, win, c1, c2, c3)

for k=1:size(input,1)
    jnk(k,:)=SmoothN(input(k, :), win);
end
matr=jnk(:,:);
inte=2.5;
y=nanmean(matr);
% y=yy(win:size(matr,2)-win);
y1=nanmean(matr)+nanstd(matr)./(size(matr, 1).^.5);
% y1=yy1(win:size(matr,2)-win);
y2=nanmean(matr)-nanstd(matr)./(size(matr, 1).^.5);
y3=y2(size(matr, 2):-1:1);
y2=y3;
Y=[y1 y2];

%x1=1:size(input,2);
%x2=size(input,2):-1:1;
x1=t;
x2=sort(t,'descend');
X=[x1 x2];
fill(X, Y, [c1/inte c2/inte c3/inte], 'LineStyle', 'none');
hold on;
%plot(1:size(input, 2), y,  'color', [c1 c2 c3],'LineWidth',1);
plot(t, y,  'color', [c1 c2 c3],'LineWidth',1);

% axis([0 size(input,2), 0 .2]);


