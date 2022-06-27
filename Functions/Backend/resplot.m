function resplot(Res,plot)


figure
for n = 1:length(plot.sigs)
    
    switch plot.sigs{n}
        
        case 'RSE'
            matplot = normdim(real(Res.Rw1w3SE(:,:,plot.tsel)));
        case 'RGSB'
            matplot = normdim(real(Res.Rw1w3GSB(:,:,plot.tsel)));
%         case 'RESA'
%             matplot = normdim((normdim(real(Res.Rw1w3ESA(:,:,plot.tsel)))));
        case 'RABS'
            if ismember('corrected',plot.opt)
                matplot = normdim(real(Res.Rw1w3ABSconv(:,:,plot.tsel))); 
            else
                matplot = normdim(real(Res.Rw1w3Absorptive(:,:,plot.tsel)));
            end
        case 'R1'
            if ismember('corrected',plot.opt)
                matplot = normdim(real(Res.Rw1w3_conv{1}(:,:,plot.tsel)));
            else
                matplot = normdim(real(Res.Rw1w3cut{1}(:,:,plot.tsel)));
            end
        case 'R2'
            if ismember('corrected',plot.opt)
                matplot = normdim(real(Res.Rw1w3_conv{2}(:,:,plot.tsel)));
            else
                matplot = normdim(real(Res.Rw1w3cut{2}(:,:,plot.tsel)));
            end
        case 'R3'
            if ismember('corrected',plot.opt)
                matplot = normdim(real(Res.Rw1w3_conv{3}(:,:,plot.tsel)));
            else
                matplot = normdim(real(Res.Rw1w3cut{3}(:,:,plot.tsel)));
            end
        case 'R4'
            if ismember('corrected',plot.opt)
                matplot = normdim(real(Res.Rw1w3_conv{4}(:,:,plot.tsel)));
            else
                matplot = normdim(real(Res.Rw1w3cut{4}(:,:,plot.tsel)));
            end
%         case 'R5'
%             matplot = normdim(real(Res.Rw1w3{5}(:,:,plot.tsel)));
%         case 'R6'
%             matplot = normdim(real(Res.Rw1w3{6}(:,:,plot.tsel)));
    end
    
    subplot(1,length(plot.sigs),n)
    if ismember('corrected',plot.opt)
        contourf(Res.w1conv./1000,Res.w3./1000,matplot,'LevelStep',0.05,'LineWidth',0.25);
    else
        contourf(Res.w1cut./1000,Res.w3cut./1000,matplot,'LevelStep',0.025,'LineWidth',0.25);
    end
    colormap(cmap2d(20));
    caxis([-1 1]);
    line('XData',[min(Res.w3cut)./1000 max(Res.w3cut./1000)],'YData',[min(Res.w3cut)./1000 max(Res.w3cut)./1000],'Linewidth',1);
    axis image
    xlim([min(Res.w1cut)./1000 max(Res.w1cut./1000)]); ylim([min(Res.w3cut)./1000 max(Res.w3cut./1000)])
    xlabel('\omega_{1} (10^{3} cm^{-1})'); ylabel('\omega_{3} (10^{3} cm^{-1})'); title(plot.sigs{n})
    set(gca, 'FontSize', 18);
end
end
