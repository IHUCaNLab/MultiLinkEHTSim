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


%%%% define some constants and puts them into a structure

c.cw_min = 15;                       % minimum contention window size
c.cw_max = (2^7)*(c.cw_min +1)-1;    % maximum contention window size
c.slot_size = 9e-6;                  % DCF slot size 
c.max_retrial = 7;    % max retrials (max number of retransmissions)

c.lower_rate = 6.5e6;     % control frames tx rate
c.sifs = 16e-6;
c.difs = 34e-6;
% N bits for service field
c.N_service = 16;
% N bits for rts/cts/ack messages
c.N_rts = 160;
c.N_cts = 112;
c.N_ack = 112;
c.N_tail = 6;

c.MAC_head = 320+c.N_service;
c.PLCP_head = 112;
% time required to transmit the phy preamble in microseconds
c.T_preamble = 40e-6;
c.T_preamble_leg = 20e-6;

c.freqs = [2.4 5 6];      % carrier frequencies
c.per = 0.1;              % packet error rate

c.ap_tx_p = 20;
c.sta_tx_p = 15;
c.guard_int = 0.8;
c.ap_ant_gain = 0;
c.sta_ant_gain = 0;
c.spatial_streams = 2;
c.bp_distance = 5;
c.n_walls = 1;
c.noise_figure =(4+10)/2 ;


c.rts = (c.T_preamble_leg + (c.N_service + c.N_rts+c.N_tail)/c.lower_rate);
c.cts = (c.T_preamble_leg + (c.N_service + c.N_cts+c.N_tail)/c.lower_rate);
c.ack = (c.T_preamble_leg + (c.N_service+c.N_ack+c.N_tail)/c.lower_rate);

c.Tc = c.rts + c.difs;

load excel_file_extended.mat

% simplifying the excel file (i.e., removing the values that are not
% necessary to speed up the mapping process.
excel_file.snr_rbir{1}(1476:end,:)=[];
excel_file.snr_rbir{2}(1641:end,:)=[];
excel_file.snr_rbir{3}(1996:end,:)=[];
excel_file.snr_rbir{4}(2314:end,:)=[];
excel_file.snr_rbir{6}(126:end,:)=[];
excel_file.snr_rbir{7}(138:end,:)=[];

c.excel_file = excel_file;

for i=1:length(excel_file.snr_per)
    c.snr_thrs(i) = interp1(excel_file.snr_per{i}(:,2), excel_file.snr_per{i}(:,1),c.per,'linear');
end

load mcs_map_80211be.mat

c.mcs_map = mcs_map;

load sinr_pdfs_rhos_-20_80db_gran_0.2__gammas_-60_150_gran_0.2.mat
rhos = -20:0.2:80;

for r=1:length(rhos)
    sinr_cdfs{r} = unique(F(:,r));
end

c.sinr_cdfs.cdf = sinr_cdfs;
c.sinr_cdfs.rhos = rhos;
c.sinr_cdfs.gammas = -60:0.2:150;


load snrs_subc_gen_model_c.mat
c.new_cdfs.cdf = vv;
c.new_cdfs.rhos = [-20:0.2:80]';

clear mcs_map rhos sinr_pdf r F sinr_cdfs excel_file vv
