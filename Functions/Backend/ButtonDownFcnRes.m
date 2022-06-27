function ButtonDownFcnRes(object, eventdata, Res, beatproc, resnum)

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

if resnum == 9
    pop = 0;
    fits = 0;
    isos = 0;
    
    for k = 1:6
        pop = pop+real(Res.Rw1w3cut{1,k});
        fits = fits+real(beatproc{1,k}.rfits);
        isos = isos+real(beatproc{1,k}.riso);
    end
    pwr = abs(beatproc{1,resnum}.dataw1w2w3);
else
pop = real(Res.Rw1w3cut{1,resnum});
pwr = abs(beatproc{1,resnum}.dataw1w2w3);
fits = (beatproc{1,resnum}.rfits);
isos = (beatproc{1,resnum}.riso);

end
pop_check(1,:) = pop(idxrow,idxcol,:);
fit_check(1,:) = fits(idxrow,idxcol,:);
residual_check(1,:) = isos(idxrow,idxcol,:);
pwr_check(1,:) = pwr(idxrow,idxcol,:);
w_beat = Res.w2;
t2 = Res.t2;

subplot(2,2,2)
scatter(t2,pop_check,2,'filled','k')
hold on
plot(t2,fit_check,'r','LineWidth',1.5)
title(c_point(1,1))
box on
xlabel('t{2} (fs)');
set(gca,'FontSize',18)
hold off

subplot(2,2,3)
plot(t2,residual_check,'k')
title(c_point(1,2))
box on
xlabel('t_{2} (fs)');
set(gca,'FontSize',18)

subplot(2,2,4)
plot(w_beat,abs(pwr_check),'k','Linewidth',2)
title(c_point(1,2))
xlim([-3500 3500])
xlabel('\omega_{beat} (cm^{-1})');
box on
set(gca, 'FontSize', 18);

pause=1; %place debugging marker here if would like to analyze workspace produced from this function

end





        