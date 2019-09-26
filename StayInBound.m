
function [ value ] = StayInBound( value )

if value >1
    value = 1;
end

if value < -1
    value = -1;
end
    
end

