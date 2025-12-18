module CondensateScaling

export mean_condensate_size

using Statistics: mean

function mean_condensate_size(lat, threshold::Real)
    return mean(lat[lat .> threshold])
end

end