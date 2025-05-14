function [SNRE] = compute_snreff(SNRs, MCS, c)
%COMPUTE_SNREFF  Returns the SNR effective given a matrix of SNRs and an
%MCS index.
%
%[SNRE] = compute_snreff(SNRs, MCS, c)
%         where:
%         - SNRs: a Nsd x Nss SNRs matrix, where Nsd is the number of
%           data subcarriers and Nss is the number of spatial streams. [dB]
%         - MCS: MCS index between 0 and 13.
%         - c: common constants structure.
%         - SNRE: SNR effective. [dB]

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


if MCS < 0 || MCS > 13
    error('Error: MCS index must be between 0 and 13.');
end

% getting the number of spatial streams and subcarriers
[Nsd,Nss] = size(SNRs);

% snr values after which the curve is stable -> no need to map
% an entry for each modulation -> 5 entries
rbir_extremes = [9.50 12.80 19.90 26.26 32 42 48];
rbir_ext_values = [1 2 4 6 8 10 12];
mod_indexes = [1 2 2 3 3 4 4 4 5 5 6 6 7 7];



rbir = zeros(Nsd*Nss,1);
SNRs_r = reshape(SNRs, [Nsd*Nss,1]);

rbir(SNRs_r >= rbir_extremes(mod_indexes(MCS+1))) = rbir_ext_values(mod_indexes(MCS+1));
mask = rbir==0;

if any(mask)
    rbir(mask)=interp1(c.excel_file.snr_rbir{mod_indexes(MCS+1)}(:, 1),c.excel_file.snr_rbir{mod_indexes(MCS+1)}(:, 2),SNRs_r(mask),'nearest', 'extrap');
end

rbirs_m = mean(rbir);

app = abs(c.excel_file.snr_rbir{mod_indexes(MCS+1)}(:,2)-rbirs_m);
[~, index] = min(app);
index=index(1);
SNRE = c.excel_file.snr_rbir{mod_indexes(MCS+1)}(index,1);

end

