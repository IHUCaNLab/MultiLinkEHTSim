function [rx_p] = get_rx_p(tx_p,x_a,x_b,y_a,y_b,z_a,z_b,ant_g_a,ant_g_b, fc, shadowing, walls, c)
%GET_RX_P  Computes the received power per antenna in dBm.
%
%[rx_p] = GET_RX_P(tx_p, x_a, x_b, y_a, y_b, z_a, z_b, ant_g_a, ant_g_b, fc, shadowing, walls, c)
%       where:
%       - tx_p     : transmission power [dBm]
%       - x_a      : transmitter x coordinate [meters]
%       - x_b      : receiver x coordinate [meters]
%       - y_a      : transmitter y coordinate [meters]
%       - y_b      : receiver y coordinate [meters]
%       - z_a      : transmitter z coordinate (height) [meters]
%       - z_b      : receiver z coordinate (height) [meters]
%       - ant_g_a  : transmitter antenna gain [dBi]
%       - ant_g_b  : receiver antenna gain [dBi]
%       - fc       : carrier frequency [GHz]
%       - shadowing: shadowing component [dB]
%       - walls    : number of walls between TX and RX
%       - c        : common constants structure
%       - rx_p     : received power per antenna [dBm]

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



    % computing the distance
    d = sqrt(abs(x_a - x_b)^2 + abs(y_a -y_b)^2 + abs(z_a - z_b)^2);

    if d > c.bp_distance
        f_part = 35 * log10(d/c.bp_distance);
    else
        f_part = 0;
    end

    % computing the path loss
    pl = 40.05 + 20*log10(fc./2.4) + 20*log10(min([c.bp_distance,d])) + f_part + 7*walls;

    tx_p_lin = 10 ^ (tx_p / 10);
    tx_p_split = tx_p_lin / c.spatial_streams ;
    tx_p_split_dbm = 10 * log10(tx_p_split);
    rx_p = tx_p_split_dbm - pl + ant_g_a + ant_g_b - shadowing;
end
