%% Global running scripts
% 

clearvars all*

zs=15;

for i = 1:5
    
    cc(:,:,i) = data_by_day(i,zs).corrcoef;
    cluster_lp(:,:,i) = data_by_day(i,zs).cluster_lp;
    cluster_d(:,:,i) = data_by_day(i,zs).cluster_d;
    cluster_r(:,:,i) = data_by_day(i,zs).cluster_r;
    cluster_u(:,:,i) = data_by_day(i,zs).cluster_u;
    all_auROC(i,:) = sort(data_by_day(day,zs).auROC);
    all_auROCd(i,:) = sort(data_by_day(day,zs).auROCd);
    all_auROCr(i,:) = sort(data_by_day(day,zs).auROCr);
    all_auROCu(i,:) = sort(data_by_day(day,zs).auROCu);
    
end

%%

figure; 

wjet = jet(11);
wjet(5,:) = [1,1,1];

subplot(2,2,1); imagesc(sum(cluster_lp,3));colormap(wjet)
title('Lever press');
subplot(2,2,2); imagesc(sum(cluster_d,3));colormap(wjet)
title('Dipper');
subplot(2,2,3); imagesc(sum(cluster_r,3));colormap(wjet)
title('Rewarded entry');
subplot(2,2,4); imagesc(sum(cluster_u,3));colormap(wjet)
title('Unrewarded entry');

figure;
subplot(1,2,1); imagesc(mean(cc,3));
xticks([1:4]); yticks([1:4]); xticklabels({'lever','dipper','rewarded','unrewarded'});
yticklabels({'lever','dipper','rewarded','unrewarded'}); 
title(sprintf('Correlation coefficients'));

subplot(1,2,2); imagesc(std(cc,[],3));
xticks([1:4]); yticks([1:4]); xticklabels({'lever','dipper','rewarded','unrewarded'});
yticklabels({'lever','dipper','rewarded','unrewarded'}); 
title(sprintf('Stdev over days for Correlation coefficients'));

set(gcf,'position',widegraph)

%
figure; 
subplot(3,2,1); plot(shiftdim(cc(1,2,:),1),'o-'); title('Lever x dipper');
subplot(3,2,2); plot(shiftdim(cc(1,3,:),1),'o-'); title('Lever x rewarded');
subplot(3,2,3); plot(shiftdim(cc(1,4,:),1),'o-'); title('Lever x unrewarded');
subplot(3,2,4); plot(shiftdim(cc(2,3,:),1),'o-'); title('Dipper x rewarded');
subplot(3,2,5); plot(shiftdim(cc(2,4,:),1),'o-'); title('Dipper x unrewarded');
subplot(3,2,6); plot(shiftdim(cc(3,4,:),1),'o-'); title('Rewarded x unrewarded');

figure; 
plot(mean(all_auROC,1),'k-'); hold on;
title('Averaged auROC for neurons (identity non-preserving)')
plot(mean(all_auROCd,1),'r-');
plot(mean(all_auROCr,1),'b-');
plot(mean(all_auROCu,1),'g-');
legend({'Lever press','Dipper','Rewarded entry','Unrewarded entry'})