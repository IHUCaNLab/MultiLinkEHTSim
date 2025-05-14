function [per] = snreff2per(snreff, mcs, c)
%SNREFF2PER  Maps specific effective SNR and MCS to PER.
%
%[PER] = SNREFF2PER(snreff, mcs, c)
%        where:
%        - snreff: effective SNR [dB]
%        - MCS   : modulation and coding scheme index (0-13)
%        - c     : common constants structure
%        - PER   : packet error rate

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


l_extremes = [-1.5 1.60 4.2 7.1 10.2 14.2 15.6 17.1 20.8 22.4 25.8 27.8 30.5 32.9];
h_extremes = [0.2 3.2   5.5 8.8 11.7 15.9 17.3 18.9 22.5 24.4 27.8 29.9 32.4 34.9];
h_values = [2.9407e-5 3.8624e-5 3.7698e-05 3.1870e-05 2.7407e-05 8.3508e-05 4.6791e-05 2.9247e-05 5.8836e-05 7.5692e-05 6.4403e-05 7.5618e-05 4.0000e-05 6.2251e-05];

per = 1;

if snreff > l_extremes(mcs+1) && snreff <= h_extremes(mcs+1)
    per=interp1(c.excel_file.snr_per{mcs+1}(:,1),c.excel_file.snr_per{mcs+1}(:,2),snreff,'nearest');
elseif snreff > h_extremes(mcs+1)
    per = h_values(mcs+1);
else
    %do nothing
end

end