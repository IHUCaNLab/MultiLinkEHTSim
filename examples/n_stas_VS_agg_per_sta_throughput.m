%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_stas_VS_agg_per_sta_throughput.m
% 
% Performs a set of simulations for networks of size 5, 10, 15 and 20, and
% for each, computes the aggregated and the per station throughput. The
% two metrics are then plotted.

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




clear; clc; close all;

addpath("../");

% load the common constants
consts;   

settings.sim_time = 2;              % simulation time to 5 seconds
settings.sta_z_range = [1 2];       % in meters
settings.ch_refresh = 100e-3;       % refresh interval for channel realizations
settings.perfect_rx = 0;            % disable perfect reception
settings.pkt_max_retx = 7;          % max packet retransmission 
settings.ch_bws = [40 40 40];       % channel bandwidths (MHz)
settings.show_scenario =0;          % disable scenario plotting
settings.ap_z = 3;                  % set AP height

settings.shad_refresh = 1;          % shadowing is refreshed every second
settings.mpdu_lifetime = inf;       % disable mpdu lifetime (set to infinity)
settings.cca_below_lims = bw2minsens(settings.ch_bws);  % CCA limits for the three bands
settings.agg_meth = 1;              % aggregation method set to TBA
settings.ap_x = 10;                 % AP x position
settings.ap_y = 10;                 % AP y position

% all packets are of size 1500 bytes, so min_pkt_size = max_pkt_size
settings.min_pkt_size = 1500*8;     
settings.max_pkt_size = 1500*8;  
% set the channel model to the statistical models for MMSE receivers
settings.model = "Stat_mod";
settings.num_agg_pkts = 50;      % max A-MPDU size set to 50 MPDUs

% repeat each simulation 5 times (not enough for convergency, but good
% enough for an example)
REALIZATIONS = 5;  
N_STAS = [5,10,15,20];
TRAFFIC_RATE = 50e6;    % each station will genereate a 50 Mbps traffic

for p = 1:length(N_STAS)
    for r = 1:REALIZATIONS

        % setting the same traffic rate to all stations
        settings.traffic_rates= ones(1,N_STAS(p)) * TRAFFIC_RATE;
        
        % generating random positions around the AP over a 10 meters radius
        % using polar coordinates
        rand_rs = 10 * sqrt(rand(N_STAS(p),1));
        rand_ths = 2*pi*rand(N_STAS(p),1);
        % and then converting them to cartesian ones
        [xes,yis] = pol2cart(rand_ths, rand_rs);
        settings.sta_exes = settings.ap_x + xes;
        settings.sta_yis = settings.ap_y + yis;
        
        %running a single simulation
        [ra, stat] = single_sim(settings,c);
        
        % computing the aggregated throughput and storing it
        agg_tp(r,p) = sum(ra.thpt);
        per_sta_tp(r,p) = mean(ra.thpt);
    end
end


% plotting of the aggregated throughput vs number of stations curve
figure(1)
subplot(2,1,1);
plot(N_STAS, mean(agg_tp)/1e6,'-*b');
xlabel("Number of stations");
ylabel("Aggregated throughput (Mbps)");
grid on;
box on;

subplot(2,1,2);
plot(N_STAS, mean(per_sta_tp)/1e6,'-*r');
xlabel("Number of stations");
ylabel("Per STA throughput (Mbps)");
grid on;
box on;



