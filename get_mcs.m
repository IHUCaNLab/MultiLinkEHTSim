function [mcs] = get_mcs(p_rx, bw)
%GET_MCS_NAIVE   Returns the MCS index given the received power and the
%channel bandwidth aiming for a PER < 10%.
%
%  [mcs] = GET_MCS(p_rx, bw)
%          where:
%          - p_rx: received power [dBm]
%          - bw  : channel bandwidth [MHz]

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
        if p_rx >= -82 && p_rx < -79
            mcs = 0;
        elseif p_rx >=-79 && p_rx <-77
            mcs = 1;
        elseif p_rx >=-77 && p_rx <-74
            mcs=2;
        elseif p_rx >=-74 && p_rx <-70
            mcs=3;
        elseif p_rx >= -70 && p_rx <-66
            mcs=4;
        elseif p_rx >= -66 && p_rx <-65
            mcs=5;
        elseif p_rx >= -65 && p_rx <-64
            mcs=6;
        elseif p_rx >= -64 && p_rx <-59
            mcs=7;
        elseif p_rx >= -59 && p_rx <-57
            mcs=8;
        elseif p_rx >= -57 && p_rx <-54
            mcs=9;
        elseif p_rx >= -54 && p_rx < -52
            mcs=10;
        elseif p_rx >= -52 && p_rx < -49
            mcs=11;
        elseif p_rx >= -49 && p_rx < -46
            mcs=12;
        elseif p_rx >= -46
            mcs=13;
        else
            error('Something went wrong!');
        end

    case 40

        if p_rx >= -79 && p_rx < -76
            mcs = 0;
        elseif p_rx >= -76 && p_rx < -74
            mcs = 1;
        elseif p_rx >= -74 && p_rx < -71
            mcs=2;
        elseif p_rx >= -71 && p_rx < -67
            mcs=3;
        elseif p_rx >= -67 && p_rx < -63
            mcs=4;
        elseif p_rx >= -63 && p_rx < -62
            mcs=5;
        elseif p_rx >= -62 && p_rx < -61
            mcs=6;
        elseif p_rx >= -61 && p_rx < -56
            mcs=7;
        elseif p_rx >= -56 && p_rx < -54
            mcs=8;
        elseif p_rx >= -54 && p_rx < -51
            mcs=9;
        elseif p_rx >= -51 && p_rx < -49
            mcs=10;
        elseif p_rx >= -49 && p_rx < -46
            mcs=11;
        elseif p_rx >= -46 && p_rx < -43
            mcs=12;
        elseif p_rx >= -43
            mcs=13;
        else
            error('Something went wrong!');
        end

    case 80
        if p_rx >= -76 && p_rx < -73
            mcs = 0;
        elseif p_rx >= -73 && p_rx < -71
            mcs = 1;
        elseif p_rx >= -71 && p_rx < -68
            mcs=2;
        elseif p_rx >= -68 && p_rx < -64
            mcs=3;
        elseif p_rx >= -64 && p_rx < -60
            mcs=4;
        elseif p_rx >= -60 && p_rx < -59
            mcs=5;
        elseif p_rx >= -59 && p_rx < -58
            mcs=6;
        elseif p_rx >= -58 && p_rx < -53
            mcs=7;
        elseif p_rx >= -53 && p_rx < -51
            mcs=8;
        elseif p_rx >= -51 && p_rx < -48
            mcs=9;
        elseif p_rx >= -48 && p_rx < -46
            mcs=10;
        elseif p_rx >= -46 && p_rx < -43
            mcs=11;
        elseif p_rx >= -43 && p_rx < -40
            mcs=12;
        elseif p_rx >= -40
            mcs=13;
        else
            error('Something went wrong!');
        end

    case 160
        if p_rx >= -73 && p_rx < -70
            mcs = 0;
        elseif p_rx >= -70 && p_rx < -68
            mcs = 1;
        elseif p_rx >= -68 && p_rx < -65
            mcs=2;
        elseif p_rx >= -65 && p_rx < -61
            mcs=3;
        elseif p_rx >= -61 && p_rx < -57
            mcs=4;
        elseif p_rx >= -57 && p_rx < -56
            mcs=5;
        elseif p_rx >= -56 && p_rx < -55
            mcs=6;
        elseif p_rx >= -55 && p_rx < -50
            mcs=7;
        elseif p_rx >= -50 && p_rx < -48
            mcs=8;
        elseif p_rx >= -48 && p_rx < -45
            mcs=9;
        elseif p_rx >= -45 && p_rx < -43
            mcs=10;
        elseif p_rx >= -43 && p_rx < -40
            mcs=11;
        elseif p_rx >= -40 && p_rx < -37
            mcs=12;
        elseif p_rx >= -37
            mcs=13;
        else
            error('Something went wrong!');
        end

    case 320
        if p_rx >= -70 && p_rx < -67
            mcs = 0;
        elseif p_rx >= -67 && p_rx < -65
            mcs = 1;
        elseif p_rx >= -65 && p_rx < -62
            mcs=2;
        elseif p_rx >= -62 && p_rx < -58
            mcs=3;
        elseif p_rx >= -58 && p_rx < -54
            mcs=4;
        elseif p_rx >= -54 && p_rx < -53
            mcs=5;
        elseif p_rx >= -53 && p_rx < -52
            mcs=6;
        elseif p_rx >= -52 && p_rx < -47
            mcs=7;
        elseif p_rx >= -47 && p_rx < -45
            mcs=8;
        elseif p_rx >= -45 && p_rx < -42
            mcs=9;
        elseif p_rx >= -42 && p_rx < -40
            mcs=10;
        elseif p_rx >= -40 && p_rx < -37
            mcs=11;
        elseif p_rx >= -37 && p_rx < -34
            mcs=12;
        elseif p_rx >= -34
            mcs=13;
        else
            error('Something went wrong!');
        end
end
end