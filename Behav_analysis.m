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
% Note: subjects 8 & 10 removed because of technical
% problems (i.e. tested with wrong parameters)


clf

Colorcode = [0.466666668653488 0.674509823322296 0.18823529779911;...
             0 0.447058826684952 0.74117648601532;...
             0.82745099067688 0.0470588244497776 0.0470588244497776];

try         

%% prepare the data    
% load raw anonymised behaviuoral data    
load('BehavPD.mat');
% remove outliers and bad trials
[nRT,nFP,~,~]=RTP(BMatFP,BMatRT);% #ok
% remove ctrl subject tested with wrong settings
nRT(:,10,:,1)=nan; %remove RTs
nFP(:,10,:,1)=nan; %remove foreperiods

%% regression
lgRT = log(nRT);%log-transform reaction times
X = bsxfun(@times,ones(size(lgRT,1),1),1:4); X = X(:);%create vector of uncertainty levels for regression
Xs = sort(repmat(X,[size(lgRT,2),1]));%same as above but collapsed across subjects


for i = 1:3; %groups
    for s = 1:size(nRT,2);%subjects
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
    meanRT = squeeze(permute((nanmean(nanmean(nRT(:,:,:,i),1),2)),[1,3,2,4]));%mean
    err=squeeze(permute(nanstd(nanmean(nRT(:,:,:,i),1),0,2),[1,3,2,4]))./sqrt(15);%SEM
  
    errorbar((1:4)',meanRT,err,'o','MarkerFaceColor','white','LineWidth',2,'MarkerSize',8,'Color',Colorcode(i,:))
end
set(gca,'FontSize',10);
set(gca,'YLim',[320,430],'YTickLabel',[],'XTickLabel',[],'XTick',(1:4),'LineWidth',2);hold on;box off;axis square;
set(gcf,'Color','white'); xlabel('uncertainty');ylabel('log(rt)'); set(gca,'FontSize',10);
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
%binnerise RT according to foreperiods
n = 7;%n-bins
mu = round(nanmean(reshape(round(nanmean(squeeze(nanmean(nanmean(nFP),2)),2)),2,2)));
sigma = round(nanmean(reshape(round(nanmean(squeeze(nanmean(nanstd(nFP),2)),2)),2,2),2));

minFP = round(min(squeeze(min(min(nFP))),[],2));
maxFP = round(max(squeeze(max(max(nFP))),[],2));

mus = [1 1 2 2]; sigs=[1 2 1 2];
for e = 1:4
    edges.ed{e} = round(norminv((0:n)/n, mu(mus(e)), sigma(sigs(e))));
    edges.ed{e}([1 end])=[minFP(e) maxFP(e)]; 
    edges.bin{e}=edges.ed{e}(1:end-1)+ round(diff(edges.ed{e})./2);
end

for group = 1:3
    for cond = 1:4
        tmpFP = nFP(:,:,cond,group);
        
      %binnerise RTs along the distribution of foreperiods  
      for bins = 1:n 
        tmpRT = lgRT(:,:,cond,group);
        idx = (tmpFP > edges.ed{cond}(bins))& (tmpFP<= edges.ed{cond}(bins+1));
        tmpFP(idx)= edges.bin{cond}(bins);
        tmpRT(~idx) = NaN;
        binMeanlgRT(bins,:,cond,group) = nanmean(tmpRT);%#ok
      end
       bFP(:,:,cond,group) = tmpFP;%#ok
       
   %Quantify the foreperiod effect as the slope of a Poisson regression
    for s = 1:14;
        tmp = binMeanlgRT(:,s,cond,group);
        b = glmfit(edges.bin{cond},exp(tmp(:)),'poisson');
        
        %store intercept & slope for statistical analysis
        intercFPEff(s,cond,group)=b(1);%#ok
        slopeFPEff(s,cond,group)=b(2);%#ok
    end  
    end
end
slopeFPEff(10,:,1) = nan;
intercFPEff(10,:,1) = nan;%#ok



%FP-EFFECT statistics
% within (reapeated): FP present for all conditions and groups (negative
% slope) Only difference between levels of uncertainty indicating that with
% higher levels of uncertainty it's hard to keep track of time. IMPORTANTLY:
% it shows that PD-OFF subjects were actively keeping track of time.

%output for SPSS
%dlmwrite('slopeFP.txt',slopeFPEff);
%dlmwrite('intercFP.txt',intercFPEff);

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

