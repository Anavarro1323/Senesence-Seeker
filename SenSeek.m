
%User Inputs
tic
for x = 1
rmin = 40;
rmax = 140; %for imfindcircles colony detection
maxpix = 450000   %for the False Positive Elimination Algoritm (FPE)
sense2 = 0.92      %For Sensitivity of imfindcircles WITH FPE - usually higher than sense1
level = 0.45;
Gens = 4

basefolder = 'C:\Users\aleja\Desktop\4B\'

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
addpath(strcat(basefolder,'Hyperlink.m'));
addpath(strcat(basefolder,'exportfig.m'));
addpath(strcat(basefolder,'crop_borders.m'));
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

%Colony Counting Wild Type
for K = [1:Gens]
   Genk = strcat(string(K),'.jpg') ;
%Image Thresholding
    for y = 1
OriginalImage = Images(K).WT ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%False Positive Elimination Algorithm
for z = 1
%FPE Setup 
for y = 1
CC1 = bwconncomp(ThreshImage); %makes a struct array of object information, saved as CC1
InLabMat = labelmatrix(CC1); %turns the CC1 structural array into a label matrix, each row is an object, each column is the label matrix is a 1D-index-of-a-pixel within the object indicated by row
RGBLabMat = label2rgb(InLabMat,'jet','w','shuffle'); %turns the sequential numbering into different colors, for visualization
end

%FPE Judgement
for y = 1
array1 = cell(CC1.NumObjects,1)    ;                                              %X represents object
for x =[1:CC1.NumObjects];                           %fills in array1 with 1 or 0, depending on whether a detected object is larger than the maxpix threshold limit. 
 if length(CC1.PixelIdxList{1,x}) < maxpix;             % compare number of pixels in object to maxpix value
 array1{x,1} = 1;                                       %Array1 is a 2D array [size is #objects = row number, 1 is column number]                                                     
 else array1{x,1} = 0;                                   %If current object has more pixels than maxpix, array1's entry for this object is 1
 end                                                        %Array1 stores whether a given object should be erased. 1 if yes, 2 if 0
 end
end

%FPE Execution
for y = 1 
 for x = [1:CC1.NumObjects];                           %If array1 is 1, delete the object, otherwise leave the object
     if array1{x,1} == 1;
     else
 InLabMat(InLabMat == x) = 0;                         %Delete the object. Within the Label matrix, a row (object) chosen to be deleted will have all of the pixel coordinates turned to 0. 
     end
 end
 FPEresult = label2rgb(InLabMat,'jet','w','shuffle');  
 imshow(FPEresult);

BWFPE =   rgb2gray(FPEresult);

ThreshImageFPE = imbinarize(BWFPE,0.95);
imshow(ThreshImageFPE);

 end
 end

 %Colony Detection
 for y = 1
[SampleC,SampleR] = imfindcircles(ThreshImageFPE,range,'ObjectPolarity','dark','Sensitivity',sense2,'Method','twostage');
figure;imshow(Images(K).WT);
viscircles(SampleC, SampleR, 'EdgeColor', 'b');

Images(K).WT_output = export_fig(gcf,strcat(OutputFolder,string(Batch1),string(Genk)),'-r500')
close all
AllQuadrants(K).WT = SampleR;

 montage({OriginalImage,ThresholdedImage,RGBLabMat,FPEresult,ThreshImageFPE,Images(K).WT_output})
 Images(K).WT_montage = export_fig(gcf,strcat(OutputFolder,string(Batch1),'montage',string(Genk)),'-r1000')
 end
end

%Colony Counting DBF4-1
for K = [1:Gens]
   Genk = strcat(string(K),'.jpg') ;
%Image Thresholding
    for y = 1
OriginalImage = Images(K).DBF4 ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%False Positive Elimination Algorithm
for z = 1
%FPE Setup 
for y = 1
CC1 = bwconncomp(ThreshImage); %makes a struct array of object information, saved as CC1
InLabMat = labelmatrix(CC1); %turns the CC1 structural array into a label matrix, each row is an object, each column is the label matrix is a 1D-index-of-a-pixel within the object indicated by row
RGBLabMat = label2rgb(InLabMat,'jet','w','shuffle'); %turns the sequential numbering into different colors, for visualization
end

%FPE Judgement
for y = 1
array1 = cell(CC1.NumObjects,1)    ;                                              %X represents object
for x =[1:CC1.NumObjects];                           %fills in array1 with 1 or 0, depending on whether a detected object is larger than the maxpix threshold limit. 
 if length(CC1.PixelIdxList{1,x}) < maxpix;             % compare number of pixels in object to maxpix value
 array1{x,1} = 1;                                       %Array1 is a 2D array [size is #objects = row number, 1 is column number]                                                     
 else array1{x,1} = 0;                                   %If current object has more pixels than maxpix, array1's entry for this object is 1
 end                                                        %Array1 stores whether a given object should be erased. 1 if yes, 2 if 0
 end
end

%FPE Execution
for y = 1 
 for x = [1:CC1.NumObjects];                           %If array1 is 1, delete the object, otherwise leave the object
     if array1{x,1} == 1;
     else
 InLabMat(InLabMat == x) = 0;                         %Delete the object. Within the Label matrix, a row (object) chosen to be deleted will have all of the pixel coordinates turned to 0. 
     end
 end
 FPEresult = label2rgb(InLabMat,'jet','w','shuffle');  
 imshow(FPEresult);

BWFPE =   rgb2gray(FPEresult);

ThreshImageFPE = imbinarize(BWFPE,0.95);
imshow(ThreshImageFPE);

 end
 end

 %Colony Detection
 for y = 1
[SampleC,SampleR] = imfindcircles(ThreshImageFPE,range,'ObjectPolarity','dark','Sensitivity',sense2,'Method','twostage');
figure;imshow(Images(K).DBF4);
viscircles(SampleC, SampleR, 'EdgeColor', 'b');

Images(K).DBF4_output = export_fig(gcf,strcat(OutputFolder,string(Batch2),string(Genk)))
close all
AllQuadrants(K).DBF4 = SampleR;

 montage({OriginalImage,ThresholdedImage,RGBLabMat,FPEresult,ThreshImageFPE,Images(K).DBF4_output})
 Images(K).DBF4_montage = export_fig(gcf,strcat(OutputFolder,string(Batch2),'montage',string(Genk)),'-r1000')
 end
end



%Colony Counting CDC13-2
for K = [1:Gens]
   Genk = strcat(string(K),'.jpg') ;
%Image Thresholding
    for y = 1
OriginalImage = Images(K).CDC13 ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%False Positive Elimination Algorithm
for z = 1
%FPE Setup 
for y = 1
CC1 = bwconncomp(ThreshImage); %makes a struct array of object information, saved as CC1
InLabMat = labelmatrix(CC1); %turns the CC1 structural array into a label matrix, each row is an object, each column is the label matrix is a 1D-index-of-a-pixel within the object indicated by row
RGBLabMat = label2rgb(InLabMat,'jet','w','shuffle'); %turns the sequential numbering into different colors, for visualization
end

%FPE Judgement
for y = 1
array1 = cell(CC1.NumObjects,1)    ;                                              %X represents object
for x =[1:CC1.NumObjects];                           %fills in array1 with 1 or 0, depending on whether a detected object is larger than the maxpix threshold limit. 
 if length(CC1.PixelIdxList{1,x}) < maxpix;             % compare number of pixels in object to maxpix value
 array1{x,1} = 1;                                       %Array1 is a 2D array [size is #objects = row number, 1 is column number]                                                     
 else array1{x,1} = 0;                                   %If current object has more pixels than maxpix, array1's entry for this object is 1
 end                                                        %Array1 stores whether a given object should be erased. 1 if yes, 2 if 0
 end
end

%FPE Execution
for y = 1 
 for x = [1:CC1.NumObjects];                           %If array1 is 1, delete the object, otherwise leave the object
     if array1{x,1} == 1;
     else
 InLabMat(InLabMat == x) = 0;                         %Delete the object. Within the Label matrix, a row (object) chosen to be deleted will have all of the pixel coordinates turned to 0. 
     end
 end
 FPEresult = label2rgb(InLabMat,'jet','w','shuffle');  
 imshow(FPEresult);

BWFPE =   rgb2gray(FPEresult);

ThreshImageFPE = imbinarize(BWFPE,0.95);
imshow(ThreshImageFPE);

 end
 end

 %Colony Detection
 for y = 1
[SampleC,SampleR] = imfindcircles(ThreshImageFPE,range,'ObjectPolarity','dark','Sensitivity',sense2,'Method','twostage');
figure;imshow(Images(K).CDC13);
viscircles(SampleC, SampleR, 'EdgeColor', 'b');

Images(K).CDC13_output = export_fig(gcf,strcat(OutputFolder,string(Batch3),string(Genk)))
close all
AllQuadrants(K).CDC13 = SampleR;


 montage({OriginalImage,ThresholdedImage,RGBLabMat,FPEresult,ThreshImageFPE,Images(K).CDC13_output})
 Images(K).CDC13_montage = export_fig(gcf,strcat(OutputFolder,string(Batch3),'montage',string(Genk)),'-r1000')
 end
end



%Colony Counting DM
for K = [1:Gens]
   Genk = strcat(string(K),'.jpg') ;
%Image Thresholding
    for y = 1
OriginalImage = Images(K).DM ;
BW =  rgb2gray(OriginalImage);
imshow(BW);
ThresholdedImage = imbinarize(BW,level);
ThreshImage = imcomplement(ThresholdedImage);
%If works, Thresholded Image should be black colonies with white background
    end



 
%False Positive Elimination Algorithm
for z = 1
%FPE Setup 
for y = 1
CC1 = bwconncomp(ThreshImage); %makes a struct array of object information, saved as CC1
InLabMat = labelmatrix(CC1); %turns the CC1 structural array into a label matrix, each row is an object, each column is the label matrix is a 1D-index-of-a-pixel within the object indicated by row
RGBLabMat = label2rgb(InLabMat,'jet','w','shuffle'); %turns the sequential numbering into different colors, for visualization
end

%FPE Judgement
for y = 1
array1 = cell(CC1.NumObjects,1)    ;                                              %X represents object
for x =[1:CC1.NumObjects];                           %fills in array1 with 1 or 0, depending on whether a detected object is larger than the maxpix threshold limit. 
 if length(CC1.PixelIdxList{1,x}) < maxpix;             % compare number of pixels in object to maxpix value
 array1{x,1} = 1;                                       %Array1 is a 2D array [size is #objects = row number, 1 is column number]                                                     
 else array1{x,1} = 0;                                   %If current object has more pixels than maxpix, array1's entry for this object is 1
 end                                                        %Array1 stores whether a given object should be erased. 1 if yes, 2 if 0
 end
end

%FPE Execution
for y = 1 
 for x = [1:CC1.NumObjects];                           %If array1 is 1, delete the object, otherwise leave the object
     if array1{x,1} == 1;
     else
 InLabMat(InLabMat == x) = 0;                         %Delete the object. Within the Label matrix, a row (object) chosen to be deleted will have all of the pixel coordinates turned to 0. 
     end
 end
 FPEresult = label2rgb(InLabMat,'jet','w','shuffle');  
 imshow(FPEresult);

BWFPE =   rgb2gray(FPEresult);

ThreshImageFPE = imbinarize(BWFPE,0.95);
imshow(ThreshImageFPE);

 end
 end

 %Colony Detection
 for y = 1
[SampleC,SampleR] = imfindcircles(ThreshImageFPE,range,'ObjectPolarity','dark','Sensitivity',sense2,'Method','twostage');
figure;imshow(Images(K).DM);
viscircles(SampleC, SampleR, 'EdgeColor', 'b');

Images(K).DM_output = export_fig(gcf,strcat(OutputFolder,string(Batch4),string(Genk)))
close all
AllQuadrants(K).DM = SampleR;

montage({OriginalImage,ThresholdedImage,RGBLabMat,FPEresult,ThreshImageFPE,Images(K).DM_output})
 Images(K).DM_montage = export_fig(gcf,strcat(OutputFolder,string(Batch4),'montage',string(Genk)),'-r1000')

 end
end

%Histogram Plot WT
for x= [1:Gens]
    f = figure;
histogram(AllQuadrants(x).WT,20,'BinLimits',[rmin rmax]);
Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
MainTitle = title({strcat(batch1 ,' Histogram Gen ',string(x))});
saveas(gcf, strcat(HistogramsFolder,batch1, 'Histogram Gen ',string(x)));
saveas(gcf, strcat(HistogramsFolder,batch1,'Histogram Gen ',string(x),'.jpg'))
close all
end


%Histogram Plot DBF4-1
for x= [1:Gens]
    f = figure;
histogram(AllQuadrants(x).DBF4,20,'BinLimits',[rmin rmax]);
Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
MainTitle = title({strcat(batch2 ,' Histogram Gen ',string(x))});
saveas(gcf, strcat(HistogramsFolder,batch2,'Histogram Gen',string(x)));
saveas(gcf, strcat(HistogramsFolder,batch2,'Histogram Gen ',string(x),'.jpg'));
close all
end


%Histogram Plot CDC13
for x= [1:Gens]
    f = figure;
histogram(AllQuadrants(x).CDC13,20,'BinLimits',[rmin rmax]);
Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
MainTitle = title({strcat(batch3,' Histogram Gen ',string(x))});
saveas(gcf, strcat(HistogramsFolder,batch3,'Histogram Gen ',string(x)));
saveas(gcf, strcat(HistogramsFolder,batch3,'Histogram Gen ',string(x),'.jpg'));

close all
end


%Histogram Plot DM
for x= [1:Gens]
    f = figure;
histogram(AllQuadrants(x).DM,20,'BinLimits',[rmin rmax]);
Yaxis = ylabel({'Colony Counts'});                      % Create ylabel
Xaxis = xlabel({'Colony Size (Radius in Pixels)'});
MainTitle = title({strcat(batch4 ,' Histogram Gen ',string(x))});
saveas(gcf, strcat(HistogramsFolder,batch4,'Histogram Gen ',string(x)));
saveas(gcf, strcat(HistogramsFolder,batch4,'Histogram Gen ',string(x),'.jpg'));

close all
end


%Boxplot Drawing
for x = 1   
%RawPlot WT:
WTData = padcat(AllQuadrants(:).WT)
WTplot = boxplot(WTData,'PlotStyle','compact')
ylim([rmin, rmax]);
A = get(get(gca,'children'),'children');
B = get(A,'tag');

ylabel({'Average Colony Size'});                      % Create ylabel
xlabel({'Generation Number)'});
title(strcat(batch1, ' Box Plot'));
saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot'));
saveas(gcf, strcat(BoxplotsFolder,batch1,'Boxplot','.jpg'));



%RawPlot DBF4-1
DBF4Data = padcat(AllQuadrants(:).DBF4)
DBF4plot = boxplot(DBF4Data,'PlotStyle','compact')
ylim([rmin, rmax]);
C = get(get(gca,'children'),'children');
D = get(C,'tag');

ylabel({'Average Colony Size'});                      % Create ylabel
xlabel({'Generation Number)'});
title(strcat(batch2, ' Box Plot'));
saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot'));
saveas(gcf, strcat(BoxplotsFolder,batch2,'Boxplot','.jpg'));

%RawPlot CDC13-2
CDC13Data = padcat(AllQuadrants(:).CDC13)
CDC13plot = boxplot(CDC13Data,'PlotStyle','compact')
ylim([rmin, rmax]);
E = get(get(gca,'children'),'children');
F = get(E,'tag');

ylabel({'Average Colony Size'});                      % Create ylabel
xlabel({'Generation Number)'});
title(strcat(batch3, ' Box Plot'));
saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot'));
saveas(gcf, strcat(BoxplotsFolder,batch3,'Boxplot','.jpg'));




%RawPlot DM
DMData = padcat(AllQuadrants(:).DM)
DMplot = boxplot(DMData,'PlotStyle','compact')
ylim([rmin, rmax]);
G = get(get(gca,'children'),'children');
H = get(G,'tag');

ylabel({'Average Colony Size'});                      % Create ylabel
xlabel({'Generation Number)'});
title(strcat(batch4, ' Box Plot'));
saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot'));
saveas(gcf, strcat(BoxplotsFolder,batch4,'Boxplot','.jpg'));
end



% for x = 1
% 
%     
% %Set Count Tables for coloring
% SizeTable = struct2table(count)
% SizeArray = table2array(SizeTable)
% 
% TestColor = struct2cell(MyColor)
% TestColor2 = reshape(TestColor,[Gens,1])
%     
% for x = 1 %Preset Colors
% Genend = plus(Gens,1)
% end
% 
%     
% for P = [1:Gens]
% %First
%    if(SizeArray(P) > 30)
%     MyColor(P).Mut =  "g"             %Green is above 50 data points (good) 
%    elseif(SizeArray(P,1)>15)
%        MyColor(P).Mut = "y"           % 20 data points is boundry between yellow (moderate)
%    else                                 % and red (poor)
%        MyColor(P).Mut = "r"
%    end   
% 
% end  %Assign New Colors, Mut  
%     
%     
%     
% SizeTable = struct2table(count)
% SizeArray = table2array(SizeTable)
% 
% 
% for P = [1:Gens]
% %First
%    if(SizeArray(P) > 30)
%     MyColor2(P).WT =  "g"             %Green is above 50 data points (good) 
%    elseif(SizeArray(P,1)>15)
%        MyColor2(P).WT = "y"           % 20 data points is boundry between yellow (moderate)
%    else                                 % and red (poor)
%        MyColor2(P).WT = "r"
%    end   
% 
% end %Assign New Colors, WT
% 
% TestColorWT = struct2cell(MyColor2)
% TestColorWT2 = reshape(TestColorWT,[Gens,1])
% end
% 
% 
% 
% Genend = plus(Gens,1)
% for x = [1:Gens]
% box = A(21-x);   % The 21st object is the first box
% set(box, 'Color', string(TestColor2(minus(Genend,x))));
% end
% 
% for x = [1:Gens]
% box = A(21-x);   % The 21st object is the first box
% set(box, 'Color', string(TestColor2(minus(Genend,x))));
% end
% end


toc
save('SenSeekWorkspace.mat');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%Statistical Analyses%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for x = 1
% %%%Student's T test%%%
% %%% Compare to WT of Same Generation
% [h1,p1,ci1,stats1] = ttest2(AllQuadrants(1).Mut,AllQuadrants(1).WT)
% [h2,p2,ci2,stats2] = ttest2(AllQuadrants(2).Mut,AllQuadrants(2).WT)
% [h3,p3,ci3,stats3] = ttest2(AllQuadrants(3).Mut,AllQuadrants(3).WT)
% [h4,p4,ci4,stats4] = ttest2(AllQuadrants(4).Mut,AllQuadrants(4).WT)
% %%% Compare to Same Mutant's Gen 1%%
% [h5,p5,ci5,stats5] = ttest2(AllQuadrants(2).Mut,AllQuadrants(1).Mut)
% [h6,p6,ci6,stats6] = ttest2(AllQuadrants(3).Mut,AllQuadrants(1).Mut)
% [h7,p7,ci7,stats7] = ttest2(AllQuadrants(4).Mut,AllQuadrants(1).Mut)
% 
% [h8,p8,ci8,stats8] = ttest2(AllQuadrants(2).WT,AllQuadrants(1).WT)
% [h9,p9,ci9,stats9] = ttest2(AllQuadrants(3).WT,AllQuadrants(1).WT)
% [h10,p10,ci10,stats10] = ttest2(AllQuadrants(4).WT,AllQuadrants(1).WT)
% 
% 
% % %%Anova
% % %Compare to WT of Same Generation
% % [p8,ci8,stats8] = anova1(AllQuadrants(1).Mut,AllQuadrants(1).WT)
% % [p9,ci9,stats9] = anova1(AllQuadrants(2).Mut,AllQuadrants(2).WT)
% % [p10,ci10,stats10] = anova1(AllQuadrants(3).Mut,AllQuadrants(3).WT)
% % [p11,ci11,stats11] = anova1(AllQuadrants(4).Mut,AllQuadrants(4).WT)
% % %%% Compare to Same Mutant's Gen 1%%
% % [p12,ci12,stats12] = anova1(AllQuadrants(2).Mut,AllQuadrants(1).Mut)
% % [p13,ci13,stats13] = anova1(AllQuadrants(3).Mut,AllQuadrants(1).Mut)
% % [p14,ci14,stats14] = anova1(AllQuadrants(4).Mut,AllQuadrants(1).Mut)
% end

% %Variance and Standard Deviation%
% for x = [1:4]
% variance(x).WT = var(AllQuadrants(x).WT)
% StandardDev(x).WT = std(AllQuadrants(x).WT)
% variance(x).Mut = var(AllQuadrants(x).Mut)
% StandardDev(x).Mut = std(AllQuadrants(x).Mut)
% end














