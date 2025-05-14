function [bw] = nsd2bw(nsd)
%NSD2BW Maps NSD subcarriers to the related bandwidth in MHz.
%
%bw = NSD2BW(nsd)
%     where:
%     - nsd: number of data subcarriers
%     - bw : channel bandwidth [MHz]

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

    if nsd == 234
        bw=20;
    elseif nsd == 486
        bw=40;
    elseif nsd == 980
        bw=80;
    elseif nsd==1960
        bw=160;
    else
        error('Something went wrong! Invalid input.');
    end

end