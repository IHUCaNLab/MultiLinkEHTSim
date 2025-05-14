function [snrs] = get_snr_cdf(rho, Nsd, Nss, c)
%GET_SNR_CDF Computes a Nsd x Nss SNR matrix given a single 1-carrier SNR.
%   [SNRS] = GET_SNR_CDF(rho, Nsd, Nss, c) 
%           where:
%           - rho: single 1-carrier SNR [dB]
%           - Nsd: number of data subcarriers
%           - Nss: number of spatial streams
%           - c: common constants structure
%
%   Given rho, Nsd and Nss, a Nsd x Nss matrix containing the SNRs for each
%   subcarrier and spatial stream is computed using precomputed lookup
%   tables. This method is much faster than the conventional MMSE model
%   involving the computation of Nsd x Nss inverse matrixes, while keeping
%   a satisfactory degree of precision.
%
%   Refence paper: "On the Distribution of SINR for MMSE MIMO Systems"
%                  doi: doi.org/10.1109/TCOMM.2019.2900678

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

edges = [-Inf, mean([c.sinr_cdfs.rhos(2:end); c.sinr_cdfs.rhos(1:end-1)]), +Inf];
rho_idx = discretize(rho, edges);

rnd_values = rand([Nsd,Nss]);

snrs = interp1(c.sinr_cdfs.cdf{rho_idx},c.sinr_cdfs.gammas(1:length(c.sinr_cdfs.cdf{rho_idx}))',rnd_values,'nearest', 'extrap');

end