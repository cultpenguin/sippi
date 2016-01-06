

function eo=iseven(input_array)
% Call: eo=iseven(input_array);
% This function outputs 1 if the length of the input array is an even
% number. If on the other hand the length of the input array is odd the
% function outputs 0.
% Knud S. Cordua, 2010

if 2*floor(length(input_array)/2)==2*length(input_array)/2
    eo=1;
else
    eo=0;
end