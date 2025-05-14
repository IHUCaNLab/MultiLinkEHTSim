function [sens] = bw2minsens(bw)
%BW2MINSENS  Returns the minimum receiver sensitivity given the channel bandwidth.
%
% sens = BW2MINSENS(bw)
%        where:
%        - bw  : channel bandwidth [MHz]
%        - sens: minimum sensitivity [dBm]

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

sens = zeros(size(bw));

for i=1:numel(bw)
    switch bw(i)
        case 20
            sens(i) = -82;

        case 40
            sens(i) = -79;

        case 80
            sens(i) = -76;

        case 160
            sens(i) = -73;

        case 320
            sens(i) = -70;

        otherwise
            error("ERROR: invalid bandwidth (accepted values: 20, 40, 80, 160, 320 MHz)");
    end
end

end