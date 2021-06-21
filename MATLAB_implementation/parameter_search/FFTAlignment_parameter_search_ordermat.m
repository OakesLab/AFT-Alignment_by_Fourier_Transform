function av_ordermat = FFTAlignment_parameter_search_order(anglemat,st)

% calculate the order parameter
ordermat = NaN(size(anglemat));
for i = st+1:size(anglemat,1)-st
    for j = st+1:size(anglemat,2)-st
        temp = anglemat(i-st:i+st,j-st:j+st);
        temp2 = repmat(anglemat(i,j),2*st+1,2*st+1);
        comp = cos(temp-temp2);
        comp = comp.*comp;
        ordermat(i,j) = 2*(nanmean(comp(:))-.5);
        
        clear temp temp2 comp
    end
end

av_ordermat = nanmedian(ordermat(:));

end