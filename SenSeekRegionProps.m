
%User Inputs
tic
for x = 1
rmin = 40;
rmax = 140; %for imfindcircles colony detection
maxpix = 450000   %for the False Positive Elimination Algoritm (FPE)
sense2 = 0.92      %For Sensitivity of imfindcircles WITH FPE - usually higher than sense1
level = 0.45;
Gens = 5

repo = 'C:\Users\aleja\Desktop\2023_Data_Chapter2\github_repo'
basefolder = 'C:\Users\aleja\Desktop\2023_Data_Chapter2\7B\'

Batch1 = 'WT_Gen'
batch1 = 'WT'
Batch2 = 'DBF4-1_Gen'
batch2 = 'Dbf4-1'
Batch3 = 'CDC13-2_Gen'
batch3 = 'Cdc13-2'
Batch4 = 'DM_Gen'
batch4 = 'DM'

end

%Folder Specification
for x = 1
se = strel('disk',1); 

MidputFolder = strcat(basefolder,'Midput\');
OutputFolder = strcat(basefolder,'Output\');
HistogramsFolder = strcat(basefolder,'Histograms\');
BoxplotsFolder = strcat(basefolder,'Boxplots\')

addpath(strcat(basefolder,'print2array.m'));
addpath(strcat(basefolder,'hyperlink.m'));
addpath(strcat(basefolder,'export_fig.m'));
addpath(strcat(basefolder,'crop_borders.m'));
addpath(repo)
range = [rmin rmax]

filePattern = fullfile(MidputFolder, '*.jpg');
theFiles = dir(filePattern);
end

%Image Reading
for x = 1
%Load Batch 1 Images
for K = [1:Gens]
Genk = strcat(string(K),'.jpg') ;
Images(K).WT = imread(fullfile(MidputFolder,strcat(string(Batch1),string(Genk))));
end

%Load Batch 2 Images
for K = [1:Gens]
Genk = strcat(string(K),'.jpg') ;
Images(K).DBF4 = imread(fullfile(MidputFolder,strcat(string(Batch2),string(Genk))));
end

%Load Batch 3 Images
for K = [1:Gens]
Genk = strcat(string(K),'.jpg') ;
Images(K).CDC13 = imread(fullfile(MidputFolder,strcat(string(Batch3),string(Genk))));
end

%Load Batch 4 Images
for K = [1:Gens]
Genk = strcat(string(K),'.jpg') ;
Images(K).DM = imread(fullfile(MidputFolder,strcat(string(Batch4),string(Genk))));
end
end

%Color Map
    for y = 1

        color_map = [1, 0, 0;   % Category 1  - Max Pixel  - Red     
                     1, 0, 1;   % Category 2 - Minimum Pixel - Green  
                     1, 1, 0;   % Category 3 - Not Round Enough - Yellow    
                     1, 0, 1;   % Category 4 - "Over Round" - Magenta  
                     0, 0, 0;]   % Category 5 - Targets - Black


        color_map2 = [1, 1, 1;   % Category 1  - Max Pixel     
             1, 1, 1;   % Category 2 - Minimum Pixel   
            1, 1, 1;   % Category 3 - Not Round Enough    
             1, 1, 1;   % Category 4 - "Over Round"   
             0, 1, 1;]   % Category 5 - Targets

    end




for K = [1:Gens] 

%% %% Quadrant 1 Colony Counting (WT)

%% WT Image Assignment for Analysis
for y = 1
 %Image Thresholding
OriginalImage = Images(K).WT ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%% False Positive Elimination Algorithm

%% FPE Setup - Region Props 

stats = regionprops(ThreshImage,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
pixelcounts = arrayfun(@(s) length(s.PixelIdxList), stats);
CircleRatings = arrayfun(@(s) s.Circularity, stats);

InLabMat2 = zeros(size(ThreshImage));
for i = 1:numel(stats)
    InLabMat2(stats(i).PixelIdxList) = i;    %generates the rainbow full image
end
RGBLabMat2 = label2rgb(InLabMat2, 'jet', 'k', 'shuffle');


%% FPE Judgement
RegionCats = cell(length(stats),1);
RegionCats(:) = {0}; % Preallocate with zeros
 for x = 1:length(pixelcounts)
  
  if pixelcounts(x) > 900000 % compare number of pixels in object to maxpix valu
    RegionCats{x,1} = 1; % 'Over MaxPix';
  
  elseif pixelcounts(x) < 1000
    RegionCats{x,1} = 2; % = 'Under MinPix';
    
   elseif CircleRatings(x) < .80
    RegionCats{x,1} = 3; %'Not Round Enough';
    
   elseif CircleRatings(x) > 1.2
     RegionCats{x,1} = 4 ;%Over Round?;
    
  else 
    RegionCats{x,1} = 5 ;
    
  end

 end



%% Region Props Reset
% for i = 1:numel(stats)
%     InLabMat2(stats(i).PixelIdxList) = i;    %Resets Current Image to no
%                                               objects deleted, to change
%                                               filter parameters
% end
%% Execute FPE
  for x = 1:length(stats)
    if RegionCats{x,1} == 1
      InLabMat2(InLabMat2 == x) = 1;
    elseif RegionCats{x,1} == 2;
      InLabMat2(InLabMat2 == x) = 2;
    elseif RegionCats{x,1} == 3;
      InLabMat2(InLabMat2 == x) = 3;
    elseif RegionCats{x,1} == 4;
      InLabMat2(InLabMat2 == x) = 4;
    else
      InLabMat2(InLabMat2 == x) = 5;
    end
end




% Use label2rgb to assign a unique color to each category
FPEresult2 = label2rgb(InLabMat2, color_map);

AllImages(K).WT = FPEresult2



FPEresult3 = label2rgb(InLabMat2, color_map2);



%FPE Post Processing
BWFPE2 =   rgb2gray(FPEresult3);
ThreshImageFPE2 = imbinarize(BWFPE2,0.95);
imshow(ThreshImageFPE2);

%NewDetections
stats2 = regionprops(ThreshImageFPE2,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
TheseDetections = struct('Centroid', [], 'MajorAxisLength', [], 'Perimeter', [],'Area', [],'PixelList',[],'Circularity',[]);

for i = 1:numel(stats2)
    TheseDetections(i).Centroid = stats2(i).Centroid;
    TheseDetections(i).MajorAxisLength = stats2(i).MajorAxisLength;
    TheseDetections(i).Perimeter = stats2(i).Perimeter;
    TheseDetections(i).Area = stats2(i).Area;
    TheseDetections(i).PixelList = stats2(i).PixelList;
    TheseDetections(i).Circularity = stats2(i).Circularity;
end

AllQuadrants(K).WT = TheseDetections












%% 




%% Quadrant 2 Colony Counting (DBF4-1)
%% DBF4 Image Assignment for Analysis
for y = 1
 %Image Thresholding
OriginalImage = Images(K).DBF4 ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%% False Positive Elimination Algorithm

%% FPE Setup - Region Props 

stats = regionprops(ThreshImage,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
pixelcounts = arrayfun(@(s) length(s.PixelIdxList), stats);
CircleRatings = arrayfun(@(s) s.Circularity, stats);

InLabMat2 = zeros(size(ThreshImage));
for i = 1:numel(stats)
    InLabMat2(stats(i).PixelIdxList) = i;    %generates the rainbow full image
end
RGBLabMat2 = label2rgb(InLabMat2, 'jet', 'k', 'shuffle');


%% FPE Judgement
RegionCats = cell(length(stats),1);
RegionCats(:) = {0}; % Preallocate with zeros
 for x = 1:length(pixelcounts)
  
  if pixelcounts(x) > 900000 % compare number of pixels in object to maxpix valu
    RegionCats{x,1} = 1; % 'Over MaxPix';
  
  elseif pixelcounts(x) < 1000
    RegionCats{x,1} = 2; % = 'Under MinPix';
    
   elseif CircleRatings(x) < .80
    RegionCats{x,1} = 3; %'Not Round Enough';
    
   elseif CircleRatings(x) > 1.2
     RegionCats{x,1} = 4 ;%Over Round?;
    
  else 
    RegionCats{x,1} = 5 ;
    
  end

 end



%% Region Props Reset
% for i = 1:numel(stats)
%     InLabMat2(stats(i).PixelIdxList) = i;    %Resets Current Image to no
%                                               objects deleted, to change
%                                               filter parameters
% end
%% Execute FPE
  for x = 1:length(stats)
    if RegionCats{x,1} == 1
      InLabMat2(InLabMat2 == x) = 1;
    elseif RegionCats{x,1} == 2;
      InLabMat2(InLabMat2 == x) = 2;
    elseif RegionCats{x,1} == 3;
      InLabMat2(InLabMat2 == x) = 3;
    elseif RegionCats{x,1} == 4;
      InLabMat2(InLabMat2 == x) = 4;
    else
      InLabMat2(InLabMat2 == x) = 5;
    end
end




% Use label2rgb to assign a unique color to each category
FPEresult2 = label2rgb(InLabMat2, color_map);

AllImages(K).DBF4 = FPEresult2

FPEresult3 = label2rgb(InLabMat2, color_map2);



%FPE Post Processing
BWFPE2 =   rgb2gray(FPEresult3);
ThreshImageFPE2 = imbinarize(BWFPE2,0.95);
% imshow(ThreshImageFPE2);

%NewDetections
stats2 = regionprops(ThreshImageFPE2,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
TheseDetections = struct('Centroid', [], 'MajorAxisLength', [], 'Perimeter', [],'Area', [],'PixelList',[],'Circularity',[]);

for i = 1:numel(stats2)
    TheseDetections(i).Centroid = stats2(i).Centroid;
    TheseDetections(i).MajorAxisLength = stats2(i).MajorAxisLength;
    TheseDetections(i).Perimeter = stats2(i).Perimeter;
    TheseDetections(i).Area = stats2(i).Area;
    TheseDetections(i).PixelList = stats2(i).PixelList;
    TheseDetections(i).Circularity = stats2(i).Circularity;
end

AllQuadrants(K).DBF4 = TheseDetections










%% 





%% Quadrant 3 Colony Counting (CDC13-2)
%% CDC13 Image Assignment for Analysis
for y = 1
 %Image Thresholding
OriginalImage = Images(K).CDC13 ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%% False Positive Elimination Algorithm

%% FPE Setup - Region Props 

stats = regionprops(ThreshImage,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
pixelcounts = arrayfun(@(s) length(s.PixelIdxList), stats);
CircleRatings = arrayfun(@(s) s.Circularity, stats);

InLabMat2 = zeros(size(ThreshImage));
for i = 1:numel(stats)
    InLabMat2(stats(i).PixelIdxList) = i;    %generates the rainbow full image
end
RGBLabMat2 = label2rgb(InLabMat2, 'jet', 'k', 'shuffle');


%% FPE Judgement

RegionCats = cell(length(stats),1);
RegionCats(:) = {0}; % Preallocate with zeros
 for x = 1:length(pixelcounts)
  
  if pixelcounts(x) > 900000 % compare number of pixels in object to maxpix valu
    RegionCats{x,1} = 1; % 'Over MaxPix';
  
  elseif pixelcounts(x) < 1000
    RegionCats{x,1} = 2; % = 'Under MinPix';
    
   elseif CircleRatings(x) < .80
    RegionCats{x,1} = 3; %'Not Round Enough';
    
   elseif CircleRatings(x) > 1.2
     RegionCats{x,1} = 4 ;%Over Round?;
    
  else 
    RegionCats{x,1} = 5 ;
    
  end

 end



%% Region Props Reset
% for i = 1:numel(stats)
%     InLabMat2(stats(i).PixelIdxList) = i;    %Resets Current Image to no
%                                               objects deleted, to change
%                                               filter parameters
% end
%% Execute FPE



 for x = 1:length(stats)
    if RegionCats{x,1} == 1
      InLabMat2(InLabMat2 == x) = 1;
    elseif RegionCats{x,1} == 2;
      InLabMat2(InLabMat2 == x) = 2;
    elseif RegionCats{x,1} == 3;
      InLabMat2(InLabMat2 == x) = 3;
    elseif RegionCats{x,1} == 4;
      InLabMat2(InLabMat2 == x) = 4;
    else
      InLabMat2(InLabMat2 == x) = 5;
    end
end




% Use label2rgb to assign a unique color to each category
FPEresult2 = label2rgb(InLabMat2, color_map);

AllImages(K).CDC13 = FPEresult2

FPEresult3 = label2rgb(InLabMat2, color_map2);



%FPE Post Processing
BWFPE2 =   rgb2gray(FPEresult3);
ThreshImageFPE2 = imbinarize(BWFPE2,0.95);
% imshow(ThreshImageFPE2);

%NewDetections
stats2 = regionprops(ThreshImageFPE2,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
TheseDetections = struct('Centroid', [], 'MajorAxisLength', [], 'Perimeter', [],'Area', [],'PixelList',[],'Circularity',[]);

for i = 1:numel(stats2)
    TheseDetections(i).Centroid = stats2(i).Centroid;
    TheseDetections(i).MajorAxisLength = stats2(i).MajorAxisLength;
    TheseDetections(i).Perimeter = stats2(i).Perimeter;
    TheseDetections(i).Area = stats2(i).Area;
    TheseDetections(i).PixelList = stats2(i).PixelList;
    TheseDetections(i).Circularity = stats2(i).Circularity;
end

AllQuadrants(K).CDC13 = TheseDetections






%% 




%% %% Quadrant 4 Colony Counting (DM)

%% DM Image Assignment for Analysis
for y = 1
 %Image Thresholding
OriginalImage = Images(K).DM ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%% False Positive Elimination Algorithm

%% FPE Setup - Region Props 

stats = regionprops(ThreshImage,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
pixelcounts = arrayfun(@(s) length(s.PixelIdxList), stats);
CircleRatings = arrayfun(@(s) s.Circularity, stats);

InLabMat2 = zeros(size(ThreshImage));
for i = 1:numel(stats)
    InLabMat2(stats(i).PixelIdxList) = i;    %generates the rainbow full image
end
RGBLabMat2 = label2rgb(InLabMat2, 'jet', 'k', 'shuffle');


%% FPE Judgement
RegionCats = cell(length(stats),1);
RegionCats(:) = {0}; % Preallocate with zeros
 for x = 1:length(pixelcounts)
  
  if pixelcounts(x) > 900000 % compare number of pixels in object to maxpix valu
    RegionCats{x,1} = 1; % 'Over MaxPix';
  
  elseif pixelcounts(x) < 1000
    RegionCats{x,1} = 2; % = 'Under MinPix';
    
   elseif CircleRatings(x) < .80
    RegionCats{x,1} = 3; %'Not Round Enough';
    
   elseif CircleRatings(x) > 1.2
     RegionCats{x,1} = 4 ;%Over Round?;
    
  else 
    RegionCats{x,1} = 5 ;
    
  end

 end



%% Region Props Reset
% for i = 1:numel(stats)
%     InLabMat2(stats(i).PixelIdxList) = i;    %Resets Current Image to no
%                                               objects deleted, to change
%                                               filter parameters
% end
%% Execute FPE


 for x = 1:length(stats)
    if RegionCats{x,1} == 1
      InLabMat2(InLabMat2 == x) = 1;
    elseif RegionCats{x,1} == 2;
      InLabMat2(InLabMat2 == x) = 2;
    elseif RegionCats{x,1} == 3;
      InLabMat2(InLabMat2 == x) = 3;
    elseif RegionCats{x,1} == 4;
      InLabMat2(InLabMat2 == x) = 4;
    else
      InLabMat2(InLabMat2 == x) = 5;
    end
end




% Use label2rgb to assign a unique color to each category
FPEresult2 = label2rgb(InLabMat2, color_map);

AllImages(K).DM = FPEresult2

FPEresult3 = label2rgb(InLabMat2, color_map2);



%FPE Post Processing
BWFPE2 =   rgb2gray(FPEresult3);
ThreshImageFPE2 = imbinarize(BWFPE2,0.95);
% imshow(ThreshImageFPE2);

%NewDetections
stats2 = regionprops(ThreshImageFPE2,'Centroid','MajorAxisLength','Perimeter','Area','PixelList','PixelIdxList','Circularity');
TheseDetections = struct('Centroid', [], 'MajorAxisLength', [], 'Perimeter', [],'Area', [],'PixelList',[],'Circularity',[]);

for i = 1:numel(stats2)
    TheseDetections(i).Centroid = stats2(i).Centroid;
    TheseDetections(i).MajorAxisLength = stats2(i).MajorAxisLength;
    TheseDetections(i).Perimeter = stats2(i).Perimeter;
    TheseDetections(i).Area = stats2(i).Area;
    TheseDetections(i).PixelList = stats2(i).PixelList;
    TheseDetections(i).Circularity = stats2(i).Circularity;
end

AllQuadrants(K).DM = TheseDetections


end
close all 

toc




% % %Histogram Plot WT
% % for x= [1:Gens]
% %     f = figure;
% % histogram(AllQuadrants(x).WT,20,'BinLimits',[rmin rmax]);
% % Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
% % Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
% % MainTitle = title({strcat(batch1 ,' Histogram Gen ',string(x))});
% % saveas(gcf, strcat(HistogramsFolder,batch1, 'Histogram Gen ',string(x)));
% % saveas(gcf, strcat(HistogramsFolder,batch1,'Histogram Gen ',string(x),'.jpg'))
% % close all
% % end
% % 
% % 
% % %Histogram Plot DBF4-1
% % for x= [1:Gens]
% %     f = figure;
% % histogram(AllQuadrants(x).DBF4,20,'BinLimits',[rmin rmax]);
% % Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
% % Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
% % MainTitle = title({strcat(batch2 ,' Histogram Gen ',string(x))});
% % saveas(gcf, strcat(HistogramsFolder,batch2,'Histogram Gen',string(x)));
% % saveas(gcf, strcat(HistogramsFolder,batch2,'Histogram Gen ',string(x),'.jpg'));
% % close all
% % end
% % 
% % 
% % %Histogram Plot CDC13
% % for x= [1:Gens]
% %     f = figure;
% % histogram(AllQuadrants(x).CDC13,20,'BinLimits',[rmin rmax]);
% % Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
% % Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
% % MainTitle = title({strcat(batch3,' Histogram Gen ',string(x))});
% % saveas(gcf, strcat(HistogramsFolder,batch3,'Histogram Gen ',string(x)));
% % saveas(gcf, strcat(HistogramsFolder,batch3,'Histogram Gen ',string(x),'.jpg'));
% % 
% % close all
% % end
% % 
% % 
% % %Histogram Plot DM
% % for x= [1:Gens]
% %     f = figure;
% % histogram(AllQuadrants(x).DM,20,'BinLimits',[rmin rmax]);
% % Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
% % Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
% % MainTitle = title({strcat(batch4 ,' Histogram Gen ',string(x))});
% % saveas(gcf, strcat(HistogramsFolder,batch4,'Histogram Gen ',string(x)));
% % saveas(gcf, strcat(HistogramsFolder,batch4,'Histogram Gen ',string(x),'.jpg'));
% % 
% % close all
% % end
% % 
% % 
% % %Raw Boxplot Drawing
% % for x = 1   
% % %RawPlot WT:
% % WTData = padcat(AllQuadrants(:).WT)
% % WTplot = boxplot(WTData,'PlotStyle','compact')
% % ylim([rmin, rmax]);
% % A = get(get(gca,'children'),'children');
% % B = get(A,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch1, 'Raw Box Plot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch1,'Raw_Boxplot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch1,'Raw_Boxplot','.jpg'));
% % 
% % 
% % 
% % %RawPlot DBF4-1
% % DBF4Data = padcat(AllQuadrants(:).DBF4)
% % DBF4plot = boxplot(DBF4Data,'PlotStyle','compact')
% % ylim([rmin, rmax]);
% % C = get(get(gca,'children'),'children');
% % D = get(C,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch2, 'Raw Box Plot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch2,'Raw_Boxplot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch2,'Raw_Boxplot','.jpg'));
% % 
% % %RawPlot CDC13-2
% % CDC13Data = padcat(AllQuadrants(:).CDC13)
% % CDC13plot = boxplot(CDC13Data,'PlotStyle','compact')
% % ylim([rmin, rmax]);
% % E = get(get(gca,'children'),'children');
% % F = get(E,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch3, 'Raw Box Plot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch3,'Raw Boxplot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch3,'Raw Boxplot','.jpg'));
% % 
% % 
% % 
% % 
% % %RawPlot DM
% % DMData = padcat(AllQuadrants(:).DM)
% % DMplot = boxplot(DMData,'PlotStyle','compact')
% % ylim([rmin, rmax]);
% % G = get(get(gca,'children'),'children');
% % H = get(G,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch4, 'Raw Box Plot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch4,'Raw_Boxplot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch4,'Raw_Boxplot','.jpg'));
% % end
% % 
% % %Corrected (Compared to same Generation WT) Boxplot
% % for y = 1
% % for x = [1:Gens]
% %          Corrected(x).WT = rdivide(AllQuadrants(x).WT,mean(AllQuadrants(x).WT))  
% %     Corrected(x).DBF4 = rdivide(AllQuadrants(x).DBF4,mean(AllQuadrants(x).WT))
% %         Corrected(x).CDC13 = rdivide(AllQuadrants(x).CDC13,mean(AllQuadrants(x).WT))
% %         Corrected(x).DM = rdivide(AllQuadrants(x).DM,mean(AllQuadrants(x).WT))
% % end
% % 
% % %Corrected Boxplot Drawing
% % for x = 1   
% % %CorrectedPlot WT:
% % WTData = padcat(Corrected(:).WT)
% % WTplot = boxplot(WTData,'PlotStyle','compact')
% % A = get(get(gca,'children'),'children');
% % B = get(A,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch1, ' Box Plot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot Corrected','.jpg'));
% % 
% % 
% % 
% % %RawPlot DBF4-1
% % DBF4Data = padcat(Corrected(:).DBF4)
% % DBF4plot = boxplot(DBF4Data,'PlotStyle','compact')
% % C = get(get(gca,'children'),'children');
% % D = get(C,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch2, ' Box Plot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot Corrected','.jpg'));
% % 
% % %RawPlot CDC13-2
% % CDC13Data = padcat(Corrected(:).CDC13)
% % CDC13plot = boxplot(CDC13Data,'PlotStyle','compact')
% % E = get(get(gca,'children'),'children');
% % F = get(E,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch3, 'Corrected Box Plot'));
% % saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot Corrected','.jpg'));
% % 
% % 
% % 
% % 
% % %RawPlot DM
% % DMData = padcat(Corrected(:).DM)
% % DMplot = boxplot(DMData,'PlotStyle','compact')
% % G = get(get(gca,'children'),'children');
% % H = get(G,'tag');
% % 
% % ylabel({'Average Colony Size'});                      % Create ylabel
% % xlabel({'Generation Number)'});
% % title(strcat(batch4, ' Box Plot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot Corrected'));
% % saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot Corrected','.jpg'));
% % end
% % end
% % 
% % %Pinned - Compared to Same Genotype Boxplot
% % for y = 1
% for x = [1:Gens]
%          Pinned(x).WT = rdivide(AllQuadrants(x).WT,mean(AllQuadrants(1).WT))  
%     Pinned(x).DBF4 = rdivide(AllQuadrants(x).DBF4,mean(AllQuadrants(1).DBF4))
%        Pinned(x).CDC13 = rdivide(AllQuadrants(x).CDC13,mean(AllQuadrants(1).CDC13))
%         Pinned(x).DM = rdivide(AllQuadrants(x).DM,mean(AllQuadrants(1).DM))
% end
% 
% %Corrected Boxplot Drawing
% for x = 1   
% %CorrectedPlot WT:
% WTData = padcat(Corrected(:).WT)
% WTplot = boxplot(WTData,'PlotStyle','compact')
% A = get(get(gca,'children'),'children');
% B = get(A,'tag');
% 
% ylabel({'Average Colony Size'});                      % Create ylabel
% xlabel({'Generation Number)'});
% title(strcat(batch1, ' Box Plot Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot_Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot_Pinned','.jpg'));
% 
% 
% 
% %RawPlot DBF4-1
% DBF4Data = padcat(Corrected(:).DBF4)
% DBF4plot = boxplot(DBF4Data,'PlotStyle','compact')
% C = get(get(gca,'children'),'children');
% D = get(C,'tag');
% 
% ylabel({'Average Colony Size'});                      % Create ylabel
% xlabel({'Generation Number)'});
% title(strcat(batch2, ' Box Plot Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot_Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot_Pinned','.jpg'));
% 
% %RawPlot CDC13-2
% CDC13Data = padcat(Pinned(:).CDC13)
% CDC13plot = boxplot(CDC13Data,'PlotStyle','compact')
% E = get(get(gca,'children'),'children');
% FPEresult2 = get(E,'tag');
% 
% ylabel({'Average Colony Size'});                      % Create ylabel
% xlabel({'Generation Number)'});
% title(strcat(batch3, ' Box Plot'));
% saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot_Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot_Pinned','.jpg'));
% 
% 
% 
% 
% %RawPlot DM
% DMData = padcat(Corrected(:).DM)
% DMplot = boxplot(DMData,'PlotStyle','compact')
% G = get(get(gca,'children'),'children');
% H = get(G,'tag');
% 
% ylabel({'Average Colony Size'});                      % Create ylabel
% xlabel({'Generation Number)'});
% title(strcat(batch4, ' Box Plot Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot_Pinned'));
% saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot_Pinned','.jpg'));
% end
% end
% 
