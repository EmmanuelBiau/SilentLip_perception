clear;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; 
cd C:\ANALYSES_DATA\analyses_realign_phase\;
load ONE_TONE_ANALYSES;load TWO_TONES_ANALYSES;

subjects = 1:29;

bin = 10;
nbins = -pi:pi/bin:pi;

% save p(Tone Onsets) across participants for Chi-squared Test;
all_x1 = ONE_TONE_ANALYSES.OneTone_phase(ONE_TONE_ANALYSES.OneTone_position==1);
all_x2 = ONE_TONE_ANALYSES.OneTone_phase(ONE_TONE_ANALYSES.OneTone_position==2);
all_x3 = TWO_TONES_ANALYSES.Tone_phase(TWO_TONES_ANALYSES.Tone_position==1);
all_x4 = TWO_TONES_ANALYSES.Tone_phase(TWO_TONES_ANALYSES.Tone_position==2);

X_ttests = [];
for s = subjects
    
    X_ttests(s,1) = s; 
    
    if ismember(s,ONE_TONE_ANALYSES.PARTICIPANT) == 1
        x1 = []; x2 = []; x3 = []; x4 = [];
        x1(1,:)  = ONE_TONE_ANALYSES.OneTone_phase(ONE_TONE_ANALYSES.OneTone_position==1 & ismember(ONE_TONE_ANALYSES.PARTICIPANT,s));
        x2(1,:)  = ONE_TONE_ANALYSES.OneTone_phase(ONE_TONE_ANALYSES.OneTone_position==2 & ismember(ONE_TONE_ANALYSES.PARTICIPANT,s));
        x3(1,:)  = TWO_TONES_ANALYSES.Tone_phase(TWO_TONES_ANALYSES.Tone_position==1 & ismember(TWO_TONES_ANALYSES.PARTICIPANT,s));
        x4(1,:)  = TWO_TONES_ANALYSES.Tone_phase(TWO_TONES_ANALYSES.Tone_position==2 & ismember(TWO_TONES_ANALYSES.PARTICIPANT,s));
        X_ttests1(s,:) = histcounts(x1,nbins);
        X_ttests2(s,:)= histcounts(x2,nbins);
        X_ttests3(s,:)= histcounts(x3,nbins);
        X_ttests4(s,:)= histcounts(x4,nbins);
    elseif ismember(s,ONE_TONE_ANALYSES.PARTICIPANT) == 0
        X_ttests1(s,:) = NaN;
        X_ttests2(s,:)= NaN;
        X_ttests3(s,:)= NaN;
        X_ttests4(s,:)= NaN;
    end

end

% remove bad subjects;
bad_subjects = [1 4 25 26 27];
X_ttests1(bad_subjects,:) = [];
X_ttests2(bad_subjects,:) = [];
X_ttests3(bad_subjects,:) = [];
X_ttests4(bad_subjects,:) = [];

% Compute sum of p(Tones)/binned phase and perform t-tests on tones'onsets distributions vs. expect distribution;
Xsquare = [] ;
Xsquare(:,1) = sum(X_ttests1); 
Xsquare(:,3) = sum(X_ttests2); 
Xsquare(:,5) = sum(X_ttests3); 
Xsquare(:,7) = sum(X_ttests4); 
Xsquare(:,2) = round(mean(Xsquare(:,1)));
Xsquare(:,4) = round(mean(Xsquare(:,3)));
Xsquare(:,6) = round(mean(Xsquare(:,5)));
Xsquare(:,8) = round(mean(Xsquare(:,7)));

Xprp1 = X_ttests1./sum(X_ttests1,2);
Xprp2 = X_ttests2./sum(X_ttests2,2);
Xprp3 = X_ttests3./sum(X_ttests3,2);
Xprp4 = X_ttests4./sum(X_ttests4,2);
mean_prp = 1/size(Xprp1,2)*ones(size(Xprp1,1),1);

p_Xprp1 = [];
for t = 1:length(Xprp1(1,:))
[~,p_Xprp(1,t),~,~] = ttest(Xprp1(:,t),mean_prp(1,1),'tail','both');
[~,p_Xprp(2,t),~,~] = ttest(Xprp2(:,t),mean_prp(:,1),'tail','both');
[~,p_Xprp(3,t),~,~] = ttest(Xprp3(:,t),mean_prp(:,1),'tail','both');
[~,p_Xprp(4,t),~,~] = ttest(Xprp4(:,t),mean_prp(:,1),'tail','both');
end

corrected_p = []; c_pvalues = [];
for crt = 1:size(p_Xprp,1)
[corrected_p(crt,:), h] = bonf_holm(p_Xprp(crt,:),0.01);
[c_pvalues(crt,:), ~, ~, ~] = fdr_BH(p_Xprp(crt,:), 0.05,false);
end

%% Plot the distribution of tones along visual phase in the two conditions;
bin = 10;
nbins = -pi:pi/bin:pi;

x_labels = [];
for p = 1:length(nbins)
    if  nbins(p) == -pi
        x_labels{1,p} = '-pi';
    elseif nbins(p) == -pi/2
        x_labels{1,p} = '-pi/2';
    elseif nbins(p) == 0
        x_labels{1,p} = '0'; 
    elseif nbins(p) == pi/2
        x_labels{1,p} = 'pi/2';
    elseif nbins(p) == pi
        x_labels{1,p} = 'pi';
    else
        x_labels{1,p} = '';
    end
end

% lot distribution of Tones in the Two tones condition;
mean_Xprp3 = [mean(Xprp3),mean(Xprp3(:,1))];
mean_Xprp4 = [mean(Xprp4),mean(Xprp4(:,1))];
errhigh3 = ([std(Xprp3),std(Xprp3(:,1))])/sqrt(size(Xprp3,1));
errlow3 = ([std(Xprp3),std(Xprp3(:,1))])/sqrt(size(Xprp3,1));
errhigh4 = ([std(Xprp4),std(Xprp4(:,1))])/sqrt(size(Xprp4,1));
errlow4 = ([std(Xprp4),std(Xprp4(:,1))])/sqrt(size(Xprp4,1));

close all; figure;
subplot(2,1,1)
b3 = bar(mean_Xprp3);
b3.BarWidth = 1;
hold on
er = errorbar(1:1:length(mean_Xprp3),mean_Xprp3,errlow3,errhigh3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
line([0 22], [mean_prp(1,1) mean_prp(1,1)],'LineWidth',2,'LineStyle',':', 'Color',[1 0 0]); 
hold off
xticks(1:1:22)
xticklabels(x_labels)
ylim([0 0.15]);
title('Two Tones Condition: First Tones')
subplot(2,1,2)
b4 = bar(mean_Xprp4);
b4.BarWidth = 1;
hold on
er = errorbar(1:1:length(mean_Xprp4),mean_Xprp4,errlow4,errhigh4);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
line([0 22], [mean_prp(1,1) mean_prp(1,1)],'LineWidth',2,'LineStyle',':', 'Color',[1 0 0]); 
hold off
xticks(1:1:22)
xticklabels(x_labels)
ylim([0 0.15]);
title('Two Tones Condition: Second Tones')

% plot distribution of Tones in the One tone condition;
mean_Xprp1 = [mean(Xprp1),mean(Xprp1(:,1))];
mean_Xprp2 = [mean(Xprp2),mean(Xprp2(:,1))];
errhigh3 = ([std(Xprp1),std(Xprp1(:,1))])/sqrt(size(Xprp1,1));
errlow3 = ([std(Xprp1),std(Xprp1(:,1))])/sqrt(size(Xprp1,1));
errhigh4 = ([std(Xprp2),std(Xprp2(:,1))])/sqrt(size(Xprp2,1));
errlow4 = ([std(Xprp2),std(Xprp2(:,1))])/sqrt(size(Xprp2,1));

figure;
subplot(2,1,1)
b3 = bar(mean_Xprp1);
b3.BarWidth = 1;
hold on
er = errorbar(1:1:length(mean_Xprp1),mean_Xprp1,errlow3,errhigh3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
line([0 22], [mean_prp(1,1) mean_prp(1,1)],'LineWidth',2,'LineStyle',':', 'Color',[1 0 0]); 
hold off
xticks(1:1:22)
xticklabels(x_labels)
ylim([0 0.15]);
title('One Tone Condition: First Tones')
subplot(2,1,2)
b4 = bar(mean_Xprp2);
b4.BarWidth = 1;
hold on
er = errorbar(1:1:length(mean_Xprp2),mean_Xprp2,errlow4,errhigh4);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
line([0 22], [mean_prp(1,1) mean_prp(1,1)],'LineWidth',2,'LineStyle',':', 'Color',[1 0 0]); 
hold off
xticks(1:1:22)
xticklabels(x_labels)
ylim([0 0.15]);
title('One Tone Condition: Second Tones')
