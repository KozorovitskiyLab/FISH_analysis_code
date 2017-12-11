function  [cent, varargout]=FastPeakFind5(d, thres, filt ,edg, type, res, fid)
%License

% Copyright (c) 2012, natan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%Analyze noisy 2D images and find peaks using local maxima (1 pixel
% resolution) or weighted centroids (sub-pixel resolution).
% The code is designed to be as fast as possible, so I kept it pretty basic.
% The code assumes that the peaks are relatively sparse, test whether there
% is too much pile up and set threshold or user defined filter accordingly.
%
% How the code works:
% In theory, each peak is a smooth point spread function (SPF), like a
% Gaussian of some size, etc. In reality, there is always noise, such as
%"salt and pepper" noise, which typically has a 1 pixel variation.
% Because the peak's PSF is assumed to be larger than 1 pixel, the "true"
% local maximum of that PSF can be obtained if we can get rid of these
% single pixel noise variations. There comes medfilt2, which is a 2D median
% filter that gets rid of "salt and pepper" noise. Next we "smooth" the
% image using conv2, so that with high probability there will be only one
% pixel in each peak that will correspond to the "true" PSF local maximum.
% The weighted centroid approach uses the same image processing, with the
% difference that it just calculated the weighted centroid of each
% connected object that was obtained following the image processing.  While
% this gives sub-pixel resolution, it can miss peaks that are very close to
% each other, and runs slightly slower. Read more about how to treat these
% cases in the relevant code commentes.
%
% Inputs:
% d     The 2D data raw image - assumes a Double\Single-precision
%       floating-point, uint8 or unit16 array. Please note that the code
%       casts the raw image to uint16 if needed.  If the image dynamic range
%       is between 0 and 1, I multiplied to fit uint16. This might not be
%       optimal for generic use, so modify according to your needs.
% thres A number between 0 and max(raw_image(:)) to remove  background
% filt  A filter matrix used to smooth the image. The filter size
%       should correspond the characteristic size of the peaks
% edg   A number>1 for skipping the first few and the last few 'edge' pixels
% res   A handle that switches between two peak finding methods:
%       1 - the local maxima method (default).
%       2 - the weighted centroid sub-pixel resolution method.
%       Note that the latter method takes ~20% more time on average.
% fid   In case the user would like to save the peak positions to a file,
%       the code assumes a "fid = fopen([filename], 'w+');" line in the
%       script that uses this function.
%
%Optional Outputs:
% cent        a 1xN vector of coordinates of peaks (x1,y1,x2,y2,...
% [cent cm]   in addition to cent, cm is a binary matrix  of size(d)
%             with 1's for peak positions. (not supported in the
%             the weighted centroid sub-pixel resolution method)
%
%Example:
%
%   p=FastPeakFind(image);
%   imagesc(image); hold on
%   plot(p(1:2:end),p(2:2:end),'r+')
%
%   Adi Natan (natan@stanford.edu)
%   Ver 1.7 , Date: Oct 10th 2013
%
%% defaults
if (nargin < 1)
    d=uint16(conv2(reshape(single( 2^14*(rand(1,1024*1024)>0.99995) ),[1024 1024]) ,fspecial('gaussian', 15,3),'same')+2^8*rand(1024));
    imagesc(d);
end

if ndims(d)>2 %I added this in case one uses imread (JPG\PNG\...).
    d=uint16(rgb2gray(d));
end

if isfloat(d) %For the case the input image is double, casting to uint16 keeps enough dynamic range while speeds up the code.
    if max(d(:))<=1
        d =  uint16( d.*2^16./(max(d(:))));
    else
        d = uint16(d);
    end
end

if (nargin < 2)
    thres = (max([min(max(d,[],1))  min(max(d,[],2))])) ;
end

if (nargin < 3)
    filt = (fspecial('gaussian', 7,1)); %if needed modify the filter according to the expected peaks sizes
end

if (nargin < 4)
    edg =3;
end

if (nargin < 6)
    res = 1;
end

if (nargin < 7)
    savefileflag = false;
else
    savefileflag = true;
end

%% Analyze image
if any(d(:))  ; %for the case of non zero raw image
    
    %d = medfilt2(d,[1,1]);
    %figure, imagesc(d)
    % apply threshold
    if isa(d,'uint8')
        d=d.*uint8(d);
    else
        d=d.*uint16(d);
    end
    %figure, imagesc(d)
    if any(d(:))   ; %for the case of the image is still non zero
        %res=1;
        % smooth image
        %d=conv2(single(d),filt,'same') ;
        
        if res ==1 % switch between local maxima and sub-pixel methods
            
              % peak find - using the local maxima approach - 1 pixel resolution
                
                % d will be noisy on the edges, and also local maxima looks
                % for nearest neighbors so edge must be at least 1. We'll skip 'edge' pixels.
                sd=size(d);
                %[x y]=find(d(edg+1:sd(1)-edg,edg+1:sd(2)-edg));
                
                % initialize outputs
                cent=[];%
                cent_map=zeros(sd);
                
                %x=x+edg-1;
                %y=y+edg-1;
                for k=2:511
                for j=2:511
                    if d(k,j)>d(k-1,j-1) &&...
                            (d(k,j)>d(k-1,j)) &&...
                            (d(k,j)>d(k-1,j+1)) &&...
                            (d(k,j)>d(k,j-1)) && ...
                            (d(k,j)>d(k,j+1)) && ...
                            (d(k,j)>d(k+1,j-1)) && ...
                            (d(k,j)>d(k+1,j)) && ...
                            (d(k,j)>d(k+1,j+1));
                         cent = [cent ;  j ; k];
                        cent_map(k,j)=cent_map(k,j)+1; % if a binary matrix output is desired
                        
                        %This part accounts for saturated regions
                    elseif (d(k,j)==255) &&...
                            cent_map(k,j)<1 &&...
                            cent_map(k-1,j-1)<1 &&...
                            cent_map(k-1,j)<1 &&...
                            cent_map(k-1,j+1)<1 &&...
                            cent_map(k,j-1)<1 &&...
                            cent_map(k+1,j-1)<1; 
                            
                        i=0;
                        
                            for x=-1:1:1
                                for y=-1:1:1
                                    if 255==d(k-1,j-1);
                                        i=i+1;
                                    end
                                end
                            end
                            if i>=2;
                               cent = [cent ;  j ; k];
                               cent_map(k,j)=cent_map(k,j)+1; % if a binary matrix output is desired 
                            end
%                                         
%                             255==d(k-1,j-1) &&...
%                             
%                             (255==d(k-1,j)) &&...
%                             
%                             (255==d(k-1,j+1)) &&...
%                                 
%                             (255==d(k,j-1)) && ...
%                             
%                             (255==d(k,j+1)) && ...
%                             (255==d(k+1,j-1)) && ...
%                                
%                             (255==d(k+1,j)) && ...
%                             (255==d(k+1,j+1));
                            
                        
                    
                     
                        
                        %All these alternatives were slower...
                        %if all(reshape( d(k,j)>=d(k-1:k+1,j-1:j+1),9,1))
                        %if  d(k,j) == max(max(d((k-1):(k+1),(j-1):(j+1))))
                        %if  d(k,j)  == max(reshape(d(k,j)  >=  d(k-1:k+1,j-1:j+1),9,1))
                        
                        
                        
                    end
                end
                end
        
                
        elseif res== 2 % find weighted centroids of processed image,  sub-pixel resolution.
                   % no edg requirement needed.
                
                % get peaks areas and centroids
                stats = regionprops(logical(d),d,'Area','WeightedCentroid');
                
                % find reliable peaks by considering only peaks with an area
                % below some limit. The weighted centroid method can be not
                % accurate if peaks are very close to one another, i.e., a
                % single peak will be detected, instead of the real number
                % of peaks. This will result in a much larger area for that
                % peak. At the moment, the code ignores that peak. If that
                % happens often consider a different threshold, or return to
                % the more robust "local maxima" method.
                % To set a proper limit, inspect your data with:
                % hist([stats.Area],min([stats.Area]):max([stats.Area]));
                % to see if the limit I used (mean+2 standard deviations)
                % is an appropriate limit for your data.
                
                
                rel_peaks_vec=[stats.Area]<=mean([stats.Area])+2*std([stats.Area]);
                cent=[stats(rel_peaks_vec).WeightedCentroid]';
                new_cent=[cent(1:2:end),cent(2:2:end)];
                new_centint=uint16(new_cent);
                cent_length=length(new_cent);
                cent_map=zeros(size(d));
                for i=1:cent_length
                    cent_map(new_centint(i,2),new_centint(i,1))=1;
                
                end
        end
    
        if savefileflag
            % previous version used dlmwrite, which can be slower than  fprinf
            %             dlmwrite([filename '.txt'],[cent],   '-append', ...
            %                 'roffset', 0,   'delimiter', '\t', 'newline', 'unix');+
            
            fprintf(fid, '%f ', cent(:));
            fprintf(fid, '\n');
            
        end
        
        
    else % in case image after threshold is all zeros
        cent=[];
        cent_map=zeros(size(d));
        if nargout>1 ;  varargout{1}=cent_map; end
        return
    end
    
else % in case raw image is all zeros (dead event)
    cent=[];
    cent_map=zeros(size(d));
    if nargout>1 ;  varargout{1}=cent_map; end
    return
end

%demo mode - no input to the function
if (nargin < 1); colormap(bone);hold on; plot(cent(1:2:end),cent(2:2:end),'rs');hold off; end

% return binary mask of centroid positions if asked for
if nargout>1 ;  varargout{1}=cent_map; end

%The code below will plot the positions of the points in relationship to
%the map
scrsz = get(0,'ScreenSize');
figure('Position', [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
imagesc(d)
hold on
plot(cent(1:2:end),cent(2:2:end),'r+')

