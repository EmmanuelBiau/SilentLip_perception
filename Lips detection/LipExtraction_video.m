%% Determine visually the window containing the speaker's mouth;
clear;close; addpath D:\Birmingham_2018\NEW_video_2020\cut_videos\;
addpath D:\Birmingham_2018\NEW_video_2020\; cd D:\Birmingham_2018\NEW_video_2020\;


STIM = 551;

    %Inspection of the area of interest: You can see the information about the video;
    VR = VideoReader(['video_', num2str(STIM), '.mp4']);     
    VR.NumberOfFrames;     
    
    %Display the first frame to localize lips zone to analyse;
    frame = 1;
    video = read(VR,[1 125]);                               
    frame_inspection = video(1:720,1:1280,:,frame); 

    %Get the coordinates of your 4 points; 
    figure; imshow(frame_inspection);

    for i = 1:4

        dcm_obj = datacursormode(1);
        set(dcm_obj,'DisplayStyle','datatip','Enable','on','UpdateFcn',@myupdatefcn)
        waitforbuttonpress

        if i == 4  
         pause
         close all
        end
    end

        
%% Determine zone containing the speaker's lips;

%Enter the coordinates of the lips zone;
lips_zone = inputdlg('Enter lips_zone','Sample', [1 50]);
values = str2num(lips_zone{1});
values(1,5) = values(1,2) - values(1,1)+1;
values(1,6) = values(1,4) - values(1,3)+1;

%configure automnatically the zone to be analysed;
Lips_window = cell(1,6);
Lips_window = array2table(Lips_window, 'Variable',{'Yb','Ya','Xb','Xa','Ylim','Xlim'});
Lips_window.Yb{1} = values(1,1);
Lips_window.Ya{1} = values(1,2);
Lips_window.Xb{1} = values(1,3);
Lips_window.Xa{1} = values(1,4);
Lips_window.Ylim{1} = values(1,5);
Lips_window.Xlim{1} = values(1,6);
clear values lips_zone

%% Detect the lips movements in the zone of interest;

sz = size(video);                                                             
Mouth = zeros(Lips_window.Ylim{1},Lips_window.Xlim{1},sz(4));              
Orig = zeros (Lips_window.Ylim{1},Lips_window.Xlim{1},3,sz(4),'uint8');
Stat = zeros(4,sz(4)); 
mask = zeros(Lips_window.Ylim{1},Lips_window.Xlim{1});

%The three parameters you can adjust to improve the lips detection in the video;

lip_contrast = 194.8; 
h = fspecial('average',[5 5]); 
se = strel('disk',3);
 
%Compute the lip detection in the three dimensions: Area, Major axis and Minor axis;
pause_fig = 'on'; 
time_pause = 0.05;

for k = 1:sz(4)

    im1 = video(Lips_window.Yb{1}:Lips_window.Ya{1},Lips_window.Xb{1}:Lips_window.Xa{1},:,k); % --> you can edit 
    im1d = double(im1);
    Orig(:,:,:,k) = im1;
    im2 = (im1d(:,:,1)*170 + im1d(:,:,2)*70 + im1d(:,:,3)*70);
    im3 = sqrt(im1d(:,:,1).^2 + im1d(:,:,2).^2 + im1d(:,:,3).^2);
    %figure; imshow(im1);

    %Assume constant area;
    co = im2./im3; 
    %Here you can see filtered frame.
    %figure; imagesc(co); 

    fi = find(co > lip_contrast); 
    [~,b] = sort(co(fi));
    mask(:) = 0;
    mask(fi(b)) = 1;
    tmp = imfilter(mask,h);
    tmp(tmp < 0.5) = 0;
    tmp(tmp > 0.5) = 1;
    tmp = imclose(tmp,se);
    B = bwconvhull(tmp,'object');
    Mouth(:,:,k) = B;   
    %figure; imagesc(B); 

    st = regionprops(B,{'Area','MajorAxisLength','MinorAxisLength'});
    size_st = length(st);

    for x = 1:size_st

        %You can edit: Please see slide #8 (point 3);
        comp(x) = st(x).MajorAxisLength; 

    end

    maxval = max(comp);
    maxind = find(comp == maxval);

    %Store the lips information in the same matrix Stat;
    Stat(1,k) = st(maxind).Area;
    Stat(2,k) = st(maxind).MajorAxisLength;
    Stat(3,k) = st(maxind).MinorAxisLength;
    Stat(4,k) = maxind; % it's just an option to see the number of field (which one is taken).

    clearvars -except k VR video sz Mouth Orig Stat h se mask STIM frame_inspection Lips_window lip_contrast pause_fig time_pause

end

%Save the Area, Minor and Major axes information in a matrix;
Stat_videos = 'D:\Birmingham_2018\NEW_video_2020\Stat_videos\';
save([Stat_videos, 'Stat_video_', num2str(STIM)], 'Stat'); 

%Check on the video the lip detection and write the video with lips contour;
mouth_video = 'D:\Birmingham_2018\NEW_video_2020\mouth_videos\'; 
VO = VideoWriter([mouth_video, 'Mouth_video_', num2str(STIM)]);
open(VO); close;

    for k = 1:length(Stat)
        imshow(squeeze(Orig(:,:,:,k))); hold;
        BD = bwboundaries(squeeze(Mouth(:,:,k)),'noholes');
        maxind = Stat(4,k);
        plot(BD{maxind}(:,2),BD{maxind}(:,1),'k','LineWidth',2)
        set(gcf, 'Position', [1237 366 673 627]);
        hold; drawnow;
        cF = getframe;
        writeVideo(VO,cF)

        switch pause_fig
            case 'on' 
            pause(time_pause);
            case 'off'
            continue
        end

        isKeyPressed = ~isempty(get(gcf,'CurrentCharacter'));
             if isKeyPressed
                 break
             end

    end

%close figure at the end;
close(gcf)

%%