function s = shutter_func(m,mC)
% Shutter function
%
% Input arguments:
% m = argument
% M = number of subtime intervals
% mC = threshold
%
% Output:
% s = 1/mC if m<=mC
%   = 0 if m>mC
% 
 s = zeros(size(m));
 idx = find(m<=mC);
 s(idx) = 1/mC;

end

