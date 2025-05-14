function [index, value] = binsearch(haystack, needle)
%BINSEARCH Finds the first element of an array that is lesser or equal than 
%a value using binary search.
%
%[index, value] = BINSEARCH(haystack, needle)
%
%Finds the first element in haystack that is lesser or equal than needle 
%using binary search. haystack must be an ordered array.
%Returns the index and the value of the first element lesser or equal
%than needle.

% ------------------------------------------------------------------------
% This source code file is part of MultiLinkEHTSim.
%
% Copyright (C) 2025 Communications and Networks Laboratory, International Hellenic University
% Authors: Daniele Medda (dmedda@ihu.gr) and Athanasios Iossifides
%
% MultiLinkEHTSim is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% MultiLinkEHTSim is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along 
% with MultiLinkEHTSim. If not, you may find one here: 
% https://www.gnu.org/licenses/gpl-3.0.html


if (needle >= haystack(1) && needle <= haystack(end))

    left = 1;
    right = length(haystack);
    index = -1;

    while left <= right
        mid = floor((left + right) / 2);
        if haystack(mid) <= needle
            index = mid;
            left = mid + 1;
        else
            right = mid - 1;
        end
    end

    if index ~= -1
        value = haystack(index);
    else
        value = [];
    end
else
    if needle < haystack(1)
        index = 0;
        value = [];
    else
        index = length(haystack);
        value = [];
    end
end