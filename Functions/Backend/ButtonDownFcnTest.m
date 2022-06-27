function ButtonDownFcnTest(object, eventdata, Res, resnum)

c_point = get (gca, 'CurrentPoint');

w3 = Res.w3cut;
w1 = Res.w1cut;

rowval = c_point(1,2);                       %value of w3 to find
colval = c_point(1,1);                       %value of w1 to find
idxrow = findind(w3,rowval);
idxcol = findind(w1,colval);

% tmprow = abs(w3-rowval); 
% [idxrow idxrow] = min(tmprow);      %index of closest value
 
% tmpcol = abs(w1-colval);
% [idxcol idxcol] = min(tmpcol);      %index of closest value

pop = real(Res.Rw1w3cut{1,resnum}); %{1,resnum});

pop_check(1,:) = normdim(pop(idxrow,idxcol,:));
t2 = Res.t2;

subplot(1,2,2)
scatter(t2,pop_check,10,'filled','k')
title(c_point(1,1))
box on
ylim([min(pop_check) max(pop_check)])
xlabel('t{2} (fs)');
set(gca,'FontSize',18)
hold off

end





        