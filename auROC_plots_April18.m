%% Looking at all animals summary (7677)
zs=15;
figure(1);

daycolor = parula(6);
set(0,'DefaultLineMarkerSize',4);

for i = 1:5
    
    auROC_tmp = data_by_day7677(i,zs).auROCu;
    auROC_tmp2 = data_by_day7677(i,zs).auROCr;
    rr5_tmp = data_by_day7677(i,zs).rr_5p;
    rr95_tmp = data_by_day7677(i,zs).rr_95p;
    ur5_tmp = data_by_day7677(i,zs).ur_5p;
    ur95_tmp = data_by_day7677(i,zs).ur_95p;
    
    figure(1); 
    subplot(2,2,1); hold on; 
    plot(sort(auROC_tmp),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Unrewarded (7677)')
    
    subplot(2,2,2); hold on; 
    plot(sort(auROC_tmp2),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (7677)')
    
    subplot(2,2,3); hold on;
    patch([1:numel(ur5_tmp),numel(ur5_tmp):-1:1,1],[ur5_tmp,ur95_tmp(end:-1:1),ur5_tmp(1)],daycolor(i,:),...
        'facealpha',0.5,'edgecolor','none')
    plot(sort(auROC_tmp),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (day 7677)')
    
    subplot(2,2,4); hold on;
    patch([1:numel(rr5_tmp),numel(rr5_tmp):-1:1,1],[rr5_tmp,rr95_tmp(end:-1:1),rr5_tmp(1)],daycolor(i,:),...
        'facealpha',0.5,'edgecolor','none')
    plot(sort(auROC_tmp2),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (day 7677)')
    lgd{i} = sprintf('Day %1.0f',i);
end

widebox = [.17,.11,1.2,.65]*1000;
set(gcf,'position',widebox)
figure(1); subplot(2,2,2); legend(lgd)
saveas(gcf,'auROCr_v_u_7677.svg');


%% Looking at all animals summary (7675)
zs=15;
figure(1);

daycolor = parula(6);
set(0,'DefaultLineMarkerSize',4);

for i = 1:5
    
    auROC_tmp = data_by_day7675(i,zs).auROCu;
    auROC_tmp2 = data_by_day7675(i,zs).auROCr;
    rr5_tmp = data_by_day7675(i,zs).rr_5p;
    rr95_tmp = data_by_day7675(i,zs).rr_95p;
    ur5_tmp = data_by_day7675(i,zs).ur_5p;
    ur95_tmp = data_by_day7675(i,zs).ur_95p;
    
    figure(1); 
    subplot(2,2,1); hold on; 
    plot(sort(auROC_tmp),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Unrewarded (7675)')
    
    subplot(2,2,2); hold on; 
    plot(sort(auROC_tmp2),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (7675)')
    
    subplot(2,2,3); hold on;
    patch([1:numel(ur5_tmp),numel(ur5_tmp):-1:1,1],[ur5_tmp,ur95_tmp(end:-1:1),ur5_tmp(1)],daycolor(i,:),...
        'facealpha',0.5,'edgecolor','none')
    plot(sort(auROC_tmp),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (7675)')
    
    subplot(2,2,4); hold on;
    patch([1:numel(rr5_tmp),numel(rr5_tmp):-1:1,1],[rr5_tmp,rr95_tmp(end:-1:1),rr5_tmp(1)],daycolor(i,:),...
        'facealpha',0.5,'edgecolor','none')
    plot(sort(auROC_tmp2),'o-','markerfacecolor','w','color',daycolor(i,:)); ylim([0,1]);
    ylim([0,1]);
    xlabel('Neurons'); ylabel('auROC'); title('Rewarded (7675)')
    lgd{i} = sprintf('Day %1.0f',i);
end

widebox = [.17,.11,1.2,.65]*1000;
set(gcf,'position',widebox)
figure(1); subplot(2,2,2); legend(lgd)
saveas(gcf,'auROCr_v_u_7675.svg');

%%
figure; 

for i = 1:5
    subplot(2,5,i); imagesc(data_by_day7677(i,15).cluster_r); title(sprintf('M7677 Day %1.0f',i));
    subplot(2,5,5+i); imagesc(data_by_day7675(i,15).cluster_r); title(sprintf('M7675 Day %1.0f',i));
end

set(gcf,'position',widebox)
saveas(gcf,'clusters_by_day7677_7675.svg');

%%


figure; 

for i = 1:5
    subplot(2,5,i); imagesc(data_by_day7677(i,15).corrcoef,[-1,1]); title(sprintf('M7677 Day %1.0f',i));
    set(gca,'XTick',[1:4],'XTickLabel',{'LP','Dip','Rew','Unrew'})
    set(gca,'YTick',[1:4],'YTickLabel',{'LP','Dip','Rew','Unrew'})
    subplot(2,5,5+i); imagesc(data_by_day7675(i,15).corrcoef,[-1,1]); title(sprintf('M7675 Day %1.0f',i));
    set(gca,'XTick',[1:4],'XTickLabel',{'LP','Dip','Rew','Unrew'})
    set(gca,'YTick',[1:4],'YTickLabel',{'LP','Dip','Rew','Unrew'})
end

subplot(2,5,5); colorbar;
subplot(2,5,10); colorbar;

set(gcf,'position',widebox)
saveas(gcf,'coeffs_by_day7677_7675.svg');

%%


