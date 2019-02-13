function Behav_analysis
warning('off','all');

rng default; % For reproducibility

%% Analysis of PD dataset:
% 3 groups - Controls, Patients-off, Patients-on

%% RT analysis: RT vs Uncertainty levels
% We convert RT into log(rt) and used them to estimate intercept and slope using
% Poisson Rgression as in Banca et al 2014.
% Within group rep-ANOVA shows differences in slope and intercepts between
% on and off patients. Strikingly similar to Tomassini et al 2016.
% Note: subject 10 removed because of technical
% problems (i.e. tested with wrong parameters)


clf

Colorcode = [0.466666668653488 0.674509823322296 0.18823529779911;...
             0 0.447058826684952 0.74117648601532;...
             0.82745099067688 0.0470588244497776 0.0470588244497776];

try         

%% prepare the data    
% load raw anonymised behavioral data    
load('BehavPD.mat');%#ok
% remove outliers and bad trials
[nRT,nFP,~,~]=RTP(BMatFP,BMatRT);%#ok
% remove ctrl subject tested with wrong settings
nRT(:,10,:,1)=nan; %remove RTs
nFP(:,10,:,1)=nan; %#ok%remove foreperiods

%% regression
lgRT = log(nRT);%log-transform reaction times
X = bsxfun(@times,ones(size(lgRT,1),1),1:4); X = X(:);%create vector of uncertainty levels for regression
Xs = sort(repmat(X,[size(lgRT,2),1]));%same as above but collapsed across subjects


for i = 1:3 %groups
    for s = 1:size(nRT,2)%subjects
        tmp = nRT(:,s,:,i);
        b = glmfit(X,tmp(:),'poisson');
        
        interc(s,i)=b(1);%#ok
        slope(s,i)=b(2);%#ok
    end
   
    slope(10,1)=nan;interc(10,1)=nan;
    
    %plot regression
    tmp = nRT(:,:,:,i);
    [b,~,~] = glmfit(Xs,tmp(:),'poisson');
    yfit = glmval(b,Xs,'log');
    plot(Xs,yfit,'LineWidth',5,'Color',Colorcode(i,:));
    hold on;
    meanRT = squeeze(permute((nanmean(nanmean(nRT(:,:,:,i),2),1)),[1,3,2,4]));%mean
    
    err=squeeze(permute(nanstd(nanmean(nRT(:,:,:,i),1),0,2),[1,3,2,4]))./sqrt(15+(i~=1));%SEM
  
    errorbar((1:4)',meanRT,err,'o','MarkerFaceColor','white','LineWidth',2,'MarkerSize',8,'Color',Colorcode(i,:))
end
set(gca,'FontSize',10);
set(gca,'XLim',[0.5,4],'YTickLabel',[],'XTickLabel',[],'XTick',(1:4),'LineWidth',2);hold on;box off;axis square;
set(gcf,'Color','white'); xlabel('uncertainty');ylabel('RT'); set(gca,'FontSize',10);
...
snapnow

%plot slope and intercept
figure
subplot(1,2,1)
smallbar(nanmean(slope(:,[1 3 2])),nanstd(slope(:,[1 3 2]))./sqrt(size(lgRT,2)-1),Colorcode([1 3 2],:))
set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'LineWidth',2);hold on;box off;axis square;
set(gca,'XTickLabel',{'Control','PdOn','PdOff'},'XTick',[1 2 3]);
xlabel('group'); ylabel('slope [au]');title('slope');

subplot(1,2,2)
smallbar(nanmean(interc(:,[1 3 2])),nanstd(interc(:,[1 3 2]))./sqrt(size(lgRT,2)-1),Colorcode([1 3 2],:))
set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YLim',[5.7,6.1],'LineWidth',2);hold on;box off;axis square;
set(gca,'XTickLabel',{'Control','PdOff','PdOn'},'XTick',[1 2 3]);
xlabel('group'); ylabel('intercept [log(ms)]');title('intercept');
snapnow

%output for SPSS
%dlmwrite('BehavSlopesPD.txt',slope);
%dlmwrite('BehavIntercPD.txt',interc);


 
%% Analyse foreperiod effects
%%
N = zeros(3,1);avtmp2 =zeros(200,3);minlength = 200;
for d = 1:3
    for s = 1:size(BMatRT,2)
        for c = 1:4
            
            tmp = [BMatRT(:,s,c,d) BMatFP(:,s,c,d)];%#ok
            
            tmp = sortrows(tmp,2);
            tmp(tmp(:,1)<0.1,:)=[];
            tmp(isnan(tmp(:,1)),:)=[];
            
            
            %smooth for plotting purposes
            tmpplot = smooth(tmp(:,2),(tmp(:,1)),0.1,'rlowess');
            %statistics on raw data
            tmp2 = tmp(:,1);
            lm = fitlm(tmp2,linspace(min(tmp(:,2)),max(tmp(:,2)),length(tmp2)));
            
            
            slope(s,c,d) = lm.Coefficients.Estimate(2);
            pslope(s,c,d) = lm.Coefficients.pValue(2);%#ok
            R2(s,c,d) = lm.Rsquared.Ordinary;%#ok
            if pslope(s,c,d)>0.05; slopeC(s,c,d) = nan;  R2C(s,c,d) = nan;%#ok
            else
                slopeC(s,c,d) = slope(s,c,d);%#ok
                R2C(s,c,d) = lm.Rsquared.Ordinary;%#ok
                avtmp2(1:length(tmpplot),d) = avtmp2(1:length(tmpplot),d)+zscore(tmpplot); N(d) = N(d)+1;%this is for plotting
                minlength = min(minlength,length(tmp2));
            end
            
            pF(s,c,d) = lm.coefTest;%#ok
            
        end
    end
end


%% Compare foreperiod effects between groups
fp_stat = nanmean(slopeC,2);
%apply Bonferroni correction
alpha = 0.05/3;
%2 sample ttest - between
[h,p_ctr_off,CI,stats]=ttest2(fp_stat(:,1,1),fp_stat(:,1,2),'alpha',alpha);%#ok
[h,p_ctr_on]=ttest2(fp_stat(:,1,1),fp_stat(:,1,3),'alpha',alpha);%#ok
%paired ttest - within
[h,p_of_fon]=ttest(fp_stat(:,1,2),fp_stat(:,1,3),'alpha',alpha);%#ok

% plot foreperiod effects
figure
for iplot = 1:3
    tmpiplot = avtmp2(1:100,iplot)./N(iplot);%mean
    tmpiplot = tmpiplot-min(tmpiplot);%make all positive
    tmpiplot = tmpiplot./max(tmpiplot);%normalize with respect to max value for plotting
plot(tmpiplot,'Color',Colorcode(iplot,:),'LineWidth',3);hold on;
end

set(gca,'YTickLabel',[0 0.5 1],'XTickLabel',[0 0.5 1],'XTick',[0 50 100],'YTick',[0 0.5 1],'LineWidth',2);hold on;axis square;ylim([-0.05 1.05]);xlim([0 100])
xlabel('Normalized Foreperiod','FontSize',20,'FontWeight','Bold');ylabel('Normalized RT','FontSize',20,'FontWeight','Bold');axis square
set(gca,'FontSize', 20); set(gcf,'Color','white')


%barplots slopes of foreperiod effects
c = slopeC(:,:,1);on = slopeC(:,:,3);off = slopeC(:,:,3);
fp_eff(1,:)=[nanmean(c(:)) nanstd(c(:))/sqrt(sum(~isnan(c(:))))];
fp_eff(2,:)=[nanmean(off(:)) nanstd(off(:))/sqrt(sum(~isnan(off(:))))];
fp_eff(3,:)=[nanmean(on(:)) nanstd(on(:))/sqrt(sum(~isnan(on(:))))];
figure
smallbar(fp_eff(:,1),fp_eff(:,2),Colorcode([1 2 3],:))
set(gca,'YTickLabel',[-3 -1.5 0],'XTickLabel',[],'YTick',[-3 -1.5 0],'XTick',[1 2 3],'LineWidth',2);hold on;box off;axis square;
set(gca,'ylim',[-3 0])
set(gca,'FontSize', 20); set(gcf,'Color','white')




catch ME
    
    display('something went wrong at %s',ME.stack)
end
end



%% Ancillary functions
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
  
%robust statistics to identify outliers
noutliers = (abs(nRT)-repmat(nanmedian(nRT),[size(nRT,1),1,1]))>(1.9*repmat(mean(nRT),[size(nRT,1),1,1]));
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

