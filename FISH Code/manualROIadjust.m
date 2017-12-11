function [totMask]=manualROIadjust(DAPIMask,DAPIMaskCleanFill)

totMask=false(512);

scrsz=get(0,'ScreenSize');


figure('Position', [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3]),imagesc(DAPIMask);

figure('Position', [scrsz(3)/3 1 scrsz(3)/1.5 scrsz(4)/1.5]),imagesc(uint8(DAPIMaskCleanFill));
h=imline(gca);

BWx=createMask(h);
while sum(BWx(:))>1
totMask=totMask |BWx;
h=imline;

BWx=createMask(h);
end
