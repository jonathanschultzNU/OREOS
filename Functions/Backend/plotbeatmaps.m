function plotbeatmaps(Res, beatproc, opt, p)

w2r = round(Res.w2);
freqinds = zeros(size(opt.mapfreqs));
for i = 1:length(opt.mapfreqs)
    freqinds(i) = findind(w2r,opt.mapfreqs(i));
end

fig = figure;
    for i=1:length(freqinds)
        if ismember('corrected',opt.misc)
        beatcontours(fig,freqinds,i,opt.limcor./1000,p,Res.w1conv_beat./1000, Res.w3cut./1000, beatproc{1,opt.signal}.dataw1w2w3norm_conv(:,:,freqinds(i)),w2r(freqinds(i)),cmap2d(20))
        else
        beatcontours(fig,freqinds,i,opt.limcor./1000,p,Res.w1cut./1000, Res.w3cut./1000, beatproc{1,opt.signal}.dataw1w2w3norm(:,:,freqinds(i)),w2r(freqinds(i)),cmap2d(20))
        end
    end

end