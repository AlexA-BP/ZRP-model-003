module CondensateScaling

export mean_condensate_size, max_condensate

using Statistics: mean

function mean_condensate_size(lat, threshold::Real)
    thresholded_lat = lat[lat .> threshold]
    if isempty(thresholded_lat)
        return 0
    end
    return mean(thresholded_lat)
end

function max_condensate(lat)
    return maximum(lat)
end

end