function [rate] = mcs2rate(mcs, ss, bw, gi, c)
% MCS2RATE  Maps a specific MCS index with a specific number of spatial
% streams, bandwidth and guard interval to the related rate.
%
% rate = MCS2RATE(mcs, ss, bw, gi, c)
%        where:
%        - mcs: MCS index in range 0 - 14
%        - ss : number of spatial streams
%        - bw : bandwidth in MHz
%        - gi : guard interval in us
%        - c  : common constants structure

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

ss_shifts = [0, 14, 28, 42, 56, 70, 84, 98];

if mcs < 0 || mcs > 13
    error("Invalid MCS index!");
end

if bw == 20
    bw_shift=0;
elseif bw == 40
    bw_shift =3;
elseif bw == 80
    bw_shift = 6;
elseif bw == 160
    bw_shift = 9;
elseif bw==320
    bw_shift = 12;
else
    error("Invalid BW!");
end

if gi == 0.8
    gi_shift =0;
elseif gi == 1.6
    gi_shift = 1;
elseif gi == 3.2
    gi_shift = 2;
else
    error("Invalid guard interval!");
end


rate = c.mcs_map(mcs+ss_shifts(ss)+1,bw_shift+gi_shift+1);

end