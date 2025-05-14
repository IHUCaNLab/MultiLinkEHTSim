function [Nsd] = bw2nsd(bw)
%BW2NSD  Maps a channel bandwidth in MHz to the corresponding NSD data
%subcarriers.
% 
%Nsd = BW2NSD(bw)
%      where:
%      - bw : channel bandwidth [MHz]
%      - Nsd: number of data subcarriers

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

switch bw
    case 20
        Nsd = 234;
    case 40
        Nsd = 468;
    case 80
        Nsd = 980;
    case 160
        Nsd = 1960;
    case 320
        Nsd = 3920;
end

end