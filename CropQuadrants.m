basefolder = 'USERINPUT'           %User Input - Folder Destination.                           % "Input,Midput, and Output"                                                                      
batch = "USERINPUT";  



InputFolder = strcat(basefolder,'Input\')   %Folder should Contain 3 subfolders-      
MidputFolder = strcat(basefolder,'Midput\')  
filePattern = fullfile(InputFolder, '*.jpg');    
theFiles = dir(filePattern);


%User Inputs
Gens = 4                                                      %Number of generations, user input here
resize = 10 

%Reading Images
for x = 1
for K = [1:Gens]
                                                      %Specify name of files to search for here
    Genk = strcat(string(K),'.jpg') ;
    InputImages(K).photos =  imresize(imread(strcat(InputFolder,batch,Genk)),resize);
end 


%InputImages().photos = imrotate(InputImages().photos, 90);
end

%Center of Plate Indication by user
for k = [1:Gens]
imshow(InputImages(k).photos)
h = datacursormode                                                      %User must select center of plate here
pause(4)
InputImages(k).centers = getCursorInfo(h)
close all
end 

%Radius of plate, determined by "imdistline"
for K = [1:Gens]
imshow(InputImages(K).photos)
h = imdistline
pause(8)
InputImages(K).radii = getDistance(h);
close all
end


%PreProcess 
for K = [1:Gens]                                             %Black out every pixel outside of the plate

% Create a circle ROI and display it over the RGB image
imRGB = InputImages(K).photos;
imshow(imRGB);
roi = images.roi.Circle(gca,'Center',[InputImages(K).centers.Position(1) InputImages(K).centers.Position(2)],'Radius',InputImages(K).radii);
% Create a binary mask from the ROI
mask = createMask(roi);

% Apply the mask to the RGB image
imMasked = imRGB;
imMasked(:,:,1) = imMasked(:,:,1) .* uint8(mask);
imMasked(:,:,2) = imMasked(:,:,2) .* uint8(mask);
imMasked(:,:,3) = imMasked(:,:,3) .* uint8(mask);

% Display the masked image
imshow(imMasked);


MaskedImages(K).photos = imMasked;

 %Loads image into workspace

[Height, Width,trash] = size(imMasked);
 
%Optional chunk: shows the crop cut positions. 
%Change where the crop is with Cx and Cy
 
Cx = InputImages(K).centers.Position(1);
Cy = InputImages(K).centers.Position(2);
xsmall = minus(Cx,InputImages(K).radii);
xlarge = plus(Cx,InputImages(K).radii);
ysmall = minus(Cy,InputImages(K).radii);
ylarge = plus(Cy,InputImages(K).radii);

Q1 = imcrop(imMasked,[xsmall ysmall InputImages(K).radii InputImages(K).radii]);
Q2 = imcrop(imMasked,[Cx ysmall InputImages(K).radii InputImages(K).radii]);
Q3 = imcrop(imMasked,[xsmall Cy InputImages(K).radii InputImages(K).radii]);
Q4 = imcrop(imMasked,[Cx Cy InputImages(K).radii InputImages(K).radii]);
    
MidputImages(K).First = Q1;                                             %Cropped Images Saved to Midput folder here
imwrite(Q1,strcat(MidputFolder,'WT_Gen',string(K),'.jpg'));

MidputImages(K).Second = Q2;
imwrite(Q2,strcat(MidputFolder,'DBF4-1_Gen',string(K),'.jpg'));


MidputImages(K).Third = Q3;
imwrite(Q3,strcat(MidputFolder,'CDC13-2_Gen',string(K),'.jpg'));


MidputImages(K).Fourth = Q4;
imwrite(Q4,strcat(MidputFolder,'DM_Gen',string(K),'.jpg'));


end



save('4CropWorkspace');
