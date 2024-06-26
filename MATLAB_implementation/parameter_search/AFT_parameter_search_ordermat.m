function av_ordermat = AFT_parameter_search_ordermat(anglemat,st)

% parameters and initialisation
ordermat = NaN(size(anglemat));

% for each neighbourhood of size 2st+1 vectors
for i = st+1:size(anglemat,1)-st
    for j = st+1:size(anglemat,2)-st
        
        % compare ref vector to neighbourhood
        temp = anglemat(i-st:i+st,j-st:j+st);
        temp2 = repmat(anglemat(i,j),2*st+1,2*st+1);
        comp = cos(temp-temp2);
        comp = comp.*comp;
        
        % remove central value (comparison ref with itself)
        idx_centre = round(size(comp,1)/2);
        comp(idx_centre, idx_centre) = NaN;
        
        % calculate order parameter
        ordermat(i,j) = 2*(mean(comp(:),'omitnan')-.5);
        
        clear temp temp2 comp
    end
end

% calculate median order parameter for image
av_ordermat = median(ordermat(:),'omitnan');

end