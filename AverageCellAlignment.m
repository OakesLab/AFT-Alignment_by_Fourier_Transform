%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code can be found at https://github.com/OakesLab/FFT_Alignment
% Routine to average the alignment vectors produced from FFTAlignment.m in
% small regions.
%
% Input:
%       imname      : string of the image name to analyze (i.e. '560.tif')
%       maskname    : string of the mask indicating the regions to analyze.
%                     Image should be black background with white regions to
%                     analyze (ie. 'cellmask.tif')
%       FFTdataname : string for a .mat file that has the FFTAlignmentData
%                     variable saved from running FFTAlignment.m (i.e.
%                     'FFTactinalignment.mat')
%
% Output:
%       CellAve  : matrix with 6 columns and number of rows equal to the
%                  number of regions
%                  column (1) = region number
%                  column (2) = mean angle of alignment in that region (0 is
%                               straight down, 90 to the right, 180 straight
%                               up). Unit is degrees
%                  column (3) = horizontal position of center of mass of region
%                  column (4) = vertical position of center of mass of region
%                  column (5) = horz component of mean alignment
%                  column (6) = vertical component of mean alignment
%       CellData : struture file containing the individual alignment
%                  vectors for each region
%
% All rights and permissions belong to:
% Patrick Oakes
% poakes@gmail.com
% 09/06/2015
%
% Citation:
% Cetera et al. Nat Commun 2014; 5:1-12
% http://www.ncbi.nlm.nih.gov/pubmed/25413675
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CellAve CellData] = AverageCellAlignment(imname,maskname,FFTdataname)
% load in the local alignment data
load(FFTdataname)
pos = FFTAlignmentData.pos;
vec = FFTAlignmentData.vec;
ang = FFTAlignmentData.anglemat(:);

% load in the cell outline mask
mask = logical(imread(maskname));
masklab = bwlabel(mask);
currentmask = zeros(size(mask));

% calculate the number of cells to analyze
Ncells = max(masklab(:));

% find the points where we calculated the local alignment
pts = sub2ind(size(mask),pos(:,2),pos(:,1));
Nvec = length(pts);

% This for loop calculates the average alignment for each cell
for index = 1:Ncells
    % initialize blank mask
    currentmask = currentmask*0;
    
    % find the current cell and the points in that cell
    template = find(masklab==index);
    currentmask(template) = 1;
    
    % determine which local alignment vectors are in that cell
    keep = zeros(Nvec,1);
    for k = 1:Nvec
        if currentmask(pts(k)) == 1
            keep(k) = 1;
        else
            keep(k) = NaN;
        end
    end
    
    if nansum(keep) == 0
        fprintf('No alignment vectors in region %d\n',index)
        % data for each cell
        CellData(index,1).ang = [];
        CellData(index,1).r = [];
        CellData(index,1).c = [];
        CellData(index,1).u = [];
        CellData(index,1).v = [];
        
        % store the data as the [cell number, mean angle, mean row coordinate, mean
        % column coordinate, sin(meanangle), cos(meanangle)]
        CellAve(index,1) = index;
        CellAve(index,2) = NaN;
        CellAve(index,3) = NaN;
        CellAve(index,4) = NaN;
        CellAve(index,5) = NaN;
        CellAve(index,6) = NaN;
    else
        
        % pull out the just the data for that cell
        cellr = keep.*pos(:,2);
        cellr(isnan(cellr))=[];
        cellc = keep.*pos(:,1);
        cellc(isnan(cellc))=[];
        cellu = keep.*vec(:,1);
        cellu(isnan(cellu))=[];
        cellv = keep.*vec(:,2);
        cellv(isnan(cellv))=[];
        
        clear cellang
        % calculate the angles for the vectors
        N = length(cellu);
        for i = 1:N
            if cellv(i) > 0
                if cellv(i) > cellu(i)
                    cellang(i,1) = 180/pi*atan(cellu(i)/cellv(i));
                else
                    cellang(i,1) = 90 - 180/pi*atan(cellv(i)/cellu(i));
                end
            else
                if abs(cellv(i)) > cellu(i)
                    cellang(i,1) = 180 - 180/pi*atan(cellu(i)/abs(cellv(i)));
                else
                    cellang(i,1) = 90 + 180/pi*atan(abs(cellv(i))/cellu(i));
                end
            end
        end
        
        % determine the average angle
        a = cellang;
        N = length(a);
        k = 0:5:180;
        for p = 1:length(k)
            for i = 1:N
                if a(i) < k(p)
                    a2(i,1) = a(i)+180;
                else
                    a2(i,1) = a(i);
                end
            end
            f(p,1) = p;
            f(p,2) = mean(a2);
            f(p,3) = std(a2);
            clear a2
        end
        smin = find(f(:,3)==min(f(:,3)),1,'first');
        meancellang = f(smin,2);
        
        % make sure the angle is between 0 and 180 degrees
        if meancellang > 180
            meancellang = meancellang - 180;
        end
        if meancellang < 0
            meancellang = meancellang + 180;
        end
        
        % data for each cell
        CellData(index,1).ang = cellang;
        CellData(index,1).r = cellr;
        CellData(index,1).c = cellc;
        CellData(index,1).u = cellu;
        CellData(index,1).v = cellv;
        
        % store the data as the [cell number, mean angle, mean row coordinate, mean
        % column coordinate, sin(meanangle), cos(meanangle)]
        CellAve(index,1) = index;
        CellAve(index,2) = meancellang;
        CellAve(index,3) = round(mean(cellc));
        CellAve(index,4) = round(mean(cellr));
        CellAve(index,5) = sin(meancellang*pi/180);
        CellAve(index,6) = cos(meancellang*pi/180);
        
        clear f p a
    end
end

% Show Original Image
figure
im = imread(imname);
imshow(im,[min(im(:)) max(im(:))*.8])
hold on
% plot the average for each cell
quiver(CellAve(:,3),CellAve(:,4),CellAve(:,5),CellAve(:,6),.25,'r')

% % Plot vectors for each cell as separate colors
% co = jet(Ncells);
% for i = 1:Ncells
%     quiver(data(i).c,data(i).r,data(i).u,data(i).v,'Color',co(i,:))
% end

figure
imshow(mask)
hold on
% plot the average for each cell
quiver(CellAve(:,3),CellAve(:,4),CellAve(:,5),CellAve(:,6),.25,'r')



end
