%% PEAK DETECTION;
clearvars;clc; cd C:\ANALYSES_DATA\Signal_Infos\;
load MI_info;
    
% detect peak of MI in each dimension of lip aperture (vertical Mina, horizontal Maja and area Area);    
for p = 1:size(mi_mina_original,1)
    
   
    %detect the 3 greater MI peaks in the 1-20Hz spectrum;
    [~,locs_area] = findpeaks(mi_area_original(p,1:end),'npeaks',3,'SortStr','descend');
    [~,locs_mina] = findpeaks(mi_mina_original(p,1:end),'npeaks',3,'SortStr','descend');
    [~,locs_maja] = findpeaks(mi_maja_original(p,1:end),'npeaks',3,'SortStr','descend');

    %--- PEAKS AREA ---%
    if locs_area(1) > 3 &&  locs_area(1) < 18        
        Peaks_MI_area(p,1) = newnum(p); 
        Peaks_MI_area(p,2:8) = locs_area(1)-3:locs_area(1)+3;
        Peaks_MI_area(p,9:15) = mi_area_original(p, Peaks_MI_area(p,2:8));
      
        elseif locs_area(1) == 3        
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 locs_area(1)-2:locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 mi_area_original(p, Peaks_MI_area(p,3:8))];

        elseif locs_area(1) == 2                 
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 0 locs_area(1)-1:locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 0 mi_area_original(p, Peaks_MI_area(p,4:8))];        

        elseif locs_area(1) == 1              
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 0 0 locs_area(1):locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 0 0 mi_area_original(p, Peaks_MI_area(p,5:8))];         

        elseif locs_area(1) == 18        
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1)+2 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:7)) 0];

        elseif locs_area(1) == 19                
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1)+1 0 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:6)) 0 0];       

        elseif locs_area(1) == 20             
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1) 0 0 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:5)) 0 0 0];      
    end

    %--- PEAKS MINA ---%
    if locs_mina(1) > 3 &&  locs_mina(1) < 18        
        Peaks_MI_mina(p,1) = newnum(p); 
        Peaks_MI_mina(p,2:8) = locs_mina(1)-3:locs_mina(1)+3;
        Peaks_MI_mina(p,9:15) = mi_mina_original(p, Peaks_MI_mina(p,2:8));

        elseif locs_mina(1) == 3        
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 locs_mina(1)-2:locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 mi_mina_original(p, Peaks_MI_mina(p,3:8))];

        elseif locs_mina(1) == 2                 
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 0 locs_mina(1)-1:locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 0 mi_mina_original(p, Peaks_MI_mina(p,4:8))];        

        elseif locs_mina(1) == 1              
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 0 0 locs_mina(1):locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 0 0 mi_mina_original(p, Peaks_MI_mina(p,5:8))];         

        elseif locs_mina(1) == 18        
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1)+2 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:7)) 0];

        elseif locs_mina(1) == 19                
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1)+1 0 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:6)) 0 0];       

        elseif locs_mina(1) == 20             
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1) 0 0 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:5)) 0 0 0];     
    end

    %--- PEAKS MAJA ---%    
    if locs_maja(1) > 3 &&  locs_maja(1) < 18        
        Peaks_MI_maja(p,1) = newnum(p); 
        Peaks_MI_maja(p,2:8) = locs_maja(1)-3:locs_maja(1)+3;
        Peaks_MI_maja(p,9:15) = mi_maja_original(p, Peaks_MI_maja(p,2:8));
      
        elseif locs_maja(1) == 3        
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 locs_maja(1)-2:locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 mi_maja_original(p, Peaks_MI_maja(p,3:8))];

        elseif locs_maja(1) == 2                 
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 0 locs_maja(1)-1:locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 0 mi_maja_original(p, Peaks_MI_maja(p,4:8))];        

        elseif locs_maja(1) == 1              
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 0 0 locs_maja(1):locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 0 0 mi_maja_original(p, Peaks_MI_maja(p,5:8))];         

        elseif locs_maja(1) == 18        
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1)+2 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:7)) 0];

        elseif locs_maja(1) == 19                
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1)+1 0 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:6)) 0 0];       

        elseif locs_maja(1) == 20             
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1) 0 0 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:5)) 0 0 0];      
    end
    
    clear locs_area locs_mina locs_maja 
    
end

%% Plot Mutual information of vertical aperture (MINA): Spectrum and realigned spectrum on frequency peak of MI;
close all;figure; 
subplot(1,2,1);
spctm = 1:20;
s  = shadedErrorBar(spctm,mi_mina_original,{@mean,@std},'lineprops','-b','patchSaturation',0.5);
s2 = shadedErrorBar(spctm,mi_mina_mismatch,{@mean,@std},'lineprops','-r','patchSaturation',0.5);
set(s.mainLine, 'LineWidth',1,'DisplayName', 'Original (90videos)');
set(s.edge,'LineWidth',1,'LineStyle','none')
s.patch.FaceColor = [0.9 0.9 0.9];
set(s2.mainLine, 'LineWidth',1,'DisplayName','Valid (66videos)');
set(s2.edge,'LineWidth',1,'LineStyle','none')
s2.patch.FaceColor = [0.9 0.9 0.9];
set(gca,'XTick',1:1:20)
xlim([1 20]);
ylim([0 0.4]);
yticks([0 0.1 0.2 0.3 0.4]);
ax2 = gca;
ax2.Box = 'off';
ax2.LineWidth = 2;
ax2.TickDir = 'out';
set(ax2 ,'Layer', 'Top')
xlabel('MI Spectrum (Hz)');
ylabel('Mutual Information (MI)');
subplot(1,2,2);
plot(Peaks_MI_mina(:,9:15)', 'Color', [0.7 0.7 0.7],'LineWidth',0.25, 'LineStyle','-');
hold on
plot(mean(Peaks_MI_mina(:,9:15)), 'Color', [0 0 0],'LineWidth',2, 'LineStyle','-');
hold off
xticks([1 2 3 4 5 6 7]);
xticklabels({'fp -3','fp -2','fp -1','FreqPeak','fp +1','fp +2','fp +3'});
xlim([0.5 7.5]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
ax1 = gca;
ax1.Box = 'off';
ax1.LineWidth = 2;
ax1.TickDir = 'out';
set(ax1 ,'Layer', 'Top')
ylabel('Mutual Information (MI)');
xlabel('Realigned MI spectrum');
set(gcf,'color','w','Position',[2067 0 1502 769]);

%%