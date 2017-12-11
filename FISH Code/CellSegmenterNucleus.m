function [Ld2]=CellSegmenterNucleus(fullMask,FillMask,FillThresholds,threshold)

%Can make it so user chooses type of close disk, other parameters
%But for now I will hard code it
structElementType='disk';
structElementSize=3;
FillpixelRemovalSize=5;
fullpixelRemovalSize=5;
fullMaskBWThreshold=0;

FillMaskBWThreshold=FillThresholds(1,1)*1/256;

%process the full Mask
fullMaskNew=imsharpen(fullMask);
structElement=strel(structElementType,structElementSize);
fullMaskClose=imclose(fullMaskNew,structElement);
fullMaskNewBlur=imgaussfilt(fullMaskClose,2);
fullMaskBW=im2bw(fullMaskNewBlur,fullMaskBWThreshold);
fullMaskBWClean=bwareaopen(fullMaskBW,fullpixelRemovalSize);


%process the DAPI Mask
FillMaskBW=im2bw(FillMask,FillMaskBWThreshold);
FillMaskClean=bwareaopen(FillMaskBW,FillpixelRemovalSize);
FillMaskCleanFill=imfill(FillMaskClean,'holes');
FillMaskCleanFill=bwareaopen(FillMaskCleanFill,FillpixelRemovalSize);

totMask=manualROIadjust(FillMask,FillMaskCleanFill,fullMaskBW);
FillMaskHolder=FillMaskCleanFill;
FillMaskCleanFill(totMask>0)=0;

bw=FillMaskCleanFill;

bw3=fullMaskBWClean;



D=-bwdist(~bw);
zorromask=imextendedmin(D,2);
figure,imshowpair(bw,zorromask,'blend')
D2=imimposemin(D,zorromask);
Ld2=watershed(D2);
bw3(Ld2 == 0)=0;

%For improved speed, can comment this out if you know it is working
%properly...
figure,imshow(bw3)
figure,imshowpair(bw3,imadjust(fullMask),'blend')
figure,imshowpair(bw3,imadjust(FillMask), 'blend')
%%The line below should not be commented out
figure,imagesc(Ld2)
