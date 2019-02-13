function HDDM_analysis()


rng default; % For reproducibility

   
Colorcode = [0.466666668653488 0.674509823322296 0.18823529779911;...
             0 0.447058826684952 0.74117648601532;...
             0.82745099067688 0.0470588244497776 0.0470588244497776];

%% prepare the data    
% load raw anonymised behaviuoral data    
load('BehavPD.mat');
% remove outliers and bad trials
[nRT,nFP,~,~]=RTP(BMatFP,BMatRT);% #ok
% remove subjects tested with wrong settings

nRT(:,10,:,1)=nan;
nFP(:,10,:,1)=nan;%#ok

cols = [0 0.8 0;0 0 0.8; 0.8 0 0];

%% extract data

filename='stats_Model1.txt';
fid = fopen(filename,'rt');
tmp = textscan(fid,'%s %f %f %d %d %d %d  %d %d', 'Headerlines',1);

groups = {'Control','PD-Off','PD-On'}; 
conds = {'SL','SH','LL','LH'};
params = {'a_subj','t_subj','v_subj','t_subj'};

for parid = 1%:length(params)    
    paridx = (~cellfun(@isempty,strfind(tmp{1},params{parid})));
for gid = 1:3;
    groupidx = paridx & (~cellfun(@isempty,strfind(tmp{1},groups{gid})));
   tmp2 =[]; 
   for cid = 1:4; 
    
       condidx = groupidx & (~cellfun(@isempty,strfind(tmp{1},conds{cid})));
       tmp2(:,cid,1) = tmp{2}(condidx);  %#ok
       tmp2(:,cid,2) = tmp{2}(condidx);  %#ok       
       
   end
   
   meanv(gid).group = tmp2(:,:,1);%#ok
   
end
end
%%
for i = 1:3;
    
    RT = squeeze(nanmean(nRT(:,:,:,i),1));
    RT(isnan(RT(:,1)),:)=[];
    Par = meanv(i).group;
    for ii = 1:size(Par,1); RT(ii,:) = zscore(RT(ii,:));Par(ii,:)=zscore(Par(ii,:));end 
    
    [rho,pval]=corr(RT(:),Par(:),'Type','Spearman');
    
    subplot(3,1,i);scatter(RT(:),Par(:),20,'MarkerEdgeColor',cols(i,:),'MarkerFaceColor',cols(i,:));
    %pf = polyfit(RT(:),Par(:),1);f1 = polyval(pf,linspace(min(RT(:)),max(RT(:)),100));
    %hold on;plot(linspace(min(RT(:)),max(RT(:)),100),f1,'-r','LineWidth',2);
    axis square; box on
    xy=get(gca,'XLim');
    hold on; line(xy,xy,'LineWidth',2,'Color','k','LineStyle','--')
    if pval<0.001
        title( sprintf('%s: %0.3f; p<0.001;','\rho',rho));
    else
        title( sprintf('%s: %0.3f; p=%0.3f;','\rho',rho,pval));
    end
    set(gcf,'color','white')
end



%% Analyse slopes


mod_slopes = []; mod_int =[];
for gid = 1:3;
    tmp_data = meanv(gid).group;
    for subi = 1:size(tmp_data,1);
   
    pf = polyfit(1:4,tmp_data(subi,:),1);
    mod_slopes(subi,gid)=pf(1);%#ok
    mod_int(subi,gid)=pf(2);%#ok
    
    end
    
end

mod_slopes(end,1)=nan;mod_int(end,1)=nan;



subplot(2,1,1)
smallbar(nanmean(mod_slopes(:,[1 3 2])),nanstd(mod_slopes(:,[1 3 2]))./sqrt(16),Colorcode([1 3 2],:))
set(gca,'YTickLabel',[],'XTickLabel',{'Ctr','PDon','PDoff'},'XTick',[],'LineWidth',2);hold on;box off;axis square;
ylim([-0.02 0.1])
title('Average slope')

subplot(2,1,2)
smallbar(nanmean(mod_int(:,[1 3 2])),nanstd(mod_int(:,[1 3 2]))./sqrt(16),Colorcode([1 3 2],:))
set(gca,'YTickLabel',[],'XTickLabel',{'Ctr','PDon','PDoff'},'XTick',[1 2 3],'LineWidth',2);hold on;box off;axis square;
ylim([0.8 1.2])
set(gcf,'Color','white')
title('Average intercept')

%output for SPSS
% a = [sort(repmat([1:3]',16,1)) mod_slopes(:)];
% dlmwrite('ModSlopesPD.txt',a);
% 
% a = [sort(repmat([1:3]',16,1)) mod_int(:)];
% dlmwrite('ModIntercPD.txt',a);

%% Plot models' DIC
filename='stats_Model%s.txt';
for idm = 1:3
fid = fopen(sprintf(filename,num2str(idm)),'rt');
tmp = textscan(fid,'%s %f %f %d %d %d %d  %d %d', 'Headerlines',1);

DIC(idm) = round(tmp{2}(end-2));%#ok

end

bar(DIC,'FaceColor',[0.5 0.5 0.5]);
set(gca,'Ylim',[-24500 -24000 ])

end




% Ancillary functions

function [nRT,nFP,anticnum,numoutliers]=RTP(FP,RT)

%%%%basic preprocessing%%%%%%%%%%
%remove bad trials
nRT=RT.*1000; nFP = FP.*1000;
antic = nRT<100;
anticnum = squeeze(sum(sum(nRT<100)));
nFP(nRT<100)=NaN;
nRT(nRT<100)=NaN;%remove too early responses

%transform RTs to approac reci-normal distribution for outlier detection    
nRT = 1./nRT;
  
%robust statistics to identify ouliers
noutliers = (abs(nRT)-repmat(nanmedian(nRT),[size(nRT,1),1,1]))>(3*repmat(mad(nRT),[size(nRT,1),1,1]));
numoutliers = squeeze(sum(sum(noutliers.*not(antic))));
nRT(noutliers) = NaN;%remove outliers (x>3std)

%convert back into RTs
nRT = 1./nRT;

end


function smallbar(data,errors,cols)

if ~exist('cols','var') 
cols =  [0 0.5 0;1 0 0;...
     0.0784313753247261 0.168627455830574 0.549019634723663;...   
        ];
end

if ~exist('xlab','var')
    
    xlab = {''};%#ok
    ylab = '';%#ok
end
    h = bar(data);
set(h,'LineWidth',1);

for i = 1:3
    hold on
    bar(i,data(i),'FaceColor', cols(i,:));
    plot([i i],[data(i)-errors(i), data(i)+errors(i)],...
         '-k','LineWidth',1);
end
end



