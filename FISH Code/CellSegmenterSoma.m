function [Ld2]=CellSegmenterSoma(fullMask,FillMask,FillThresholds,threshold)

%Can make it so user chooses type of close disk, other parameters
%But for now I will hard code it
structElementType='disk';
structElementSize=3;
fillPixelRemovalSize=400;
fullpixelRemovalSize=5;
fullMaskBWThreshold=0;

fillMaskBWThreshold=FillThresholds(1,1)*1.5/256;

%process the fill Mask
fillMaskNew=imsharpen(FillMask);
structElement=strel(structElementType,structElementSize);
fillMaskClose=imclose(fillMaskNew,structElement);
fillMaskBW=im2bw(fillMaskClose,fillMaskBWThreshold);
fillMaskBWClean=bwareaopen(fillMaskBW,fillPixelRemovalSize);

totMask=manualROIadjust(FillMask,fillMaskBWClean);

fillMaskBWClean(totMask>0)=0;

[B,L]=bwboundaries(fillMaskBWClean,4,'noholes');

bw3=fillMaskBWClean;
Ld2=L;

%For improved speed, can comment this out if you know it is working
%properly...
figure,imshow(bw3)
figure,imshowpair(bw3,imadjust(fullMask),'blend')
figure,imshowpair(bw3,imadjust(FillMask), 'blend')
%The line below should not be commented out
figure,imagesc(Ld2)
