# Senesence-Seeker
This code allows cropping and a senesence seeker as of 03/25/2023

The folder containing all scripts is the base folder. 
User needs to create these folders inside of the base Folder
'Input', 'Midput', 'Output', 'Boxplots', 'Histograms'

Place all raw images of a complete petri dish (with 4 quadrants per image) into the input folder
Make sure the generation number (or whatever other sequential number you want code to iterate over) is the last character of the name of each image.

For example:

ReplicateAGen1.jpg        - A good name
Gen1ReplicateA.jpg        - NOT a good name

Batch:
"Batch refers to additional text found in the name of the file of each input image, excluding the number that will be iterated over"
for example, The "batch" for the example above is "ReplicateAGen" 

Users will first open the QuadrantCrop.m script, and update the two user inputs, explained above.


While the quadrantcrop script is running, each image will be loaded and displayed on screen. Users must click on the center of the petri dish.

Then, the images are displayed a second time, and users must place the imdistline from the center of the petri dish to the edge of the petri dish.
the length being recorded is the radius of the petri dish.


The results of the quadrantcrop script are deposited into the midput folder


Users can now run the SenSeek script.
