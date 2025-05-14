%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_known_results.m
% 
% Performs tests for the main script and save the results.
%

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

clear;
clc;
close all;

addpath("../");


stas = [5, 10, 20];
aggs = [25, 50];
trafs = [50e6, 100e6];


consts;

settings.sim_time = 1;
settings.sta_z_range = [1 2];       % in meters
settings.ch_refresh = 100e-3;       % refresh interval for channel realizations
settings.perfect_rx = 0;
settings.pkt_max_retx = 7;
settings.ch_bws = [40 40 40];       % channel bandwidths (MHz)
settings.show_scenario =0;
settings.ap_z = 3;

settings.shad_refresh = 1;
settings.mpdu_lifetime = inf;
settings.cca_below_lims = bw2minsens(settings.ch_bws);
settings.FULL_RESULTS = 0;
settings.agg_meth = 1;
settings.ap_x = 10;
settings.ap_y = 10;
settings.max_pkt_size = 1500*8;
settings.min_pkt_size = 1500*8;
settings.use_traces = 0;


n_passed=0;

l=1;
for t=1:length(trafs)
    for a=1:length(aggs)
        for s=1:length(stas)

            rng(4707);

            nstas = stas(s);
            settings.traffic_rates= ones(1,nstas) * trafs(t);
            settings.num_agg_pkts = aggs(a);
            rand_rs = 8 * sqrt(rand(nstas,1));
            rand_ths = 2*pi*rand(nstas,1);
            [xes,yis] = pol2cart(rand_ths, rand_rs);

            settings.pos_stas=[];
            settings.sta_exes = settings.ap_x + xes;
            settings.sta_yis = settings.ap_y + yis;


            [ra, stat] = single_sim(settings,c);

            dly = [];
            for k=1:nstas
                dly=[dly; mean(ra.dly_per_sta{k})];
            end
            dly=round(dly,4);
            known_results(l).thpt = ra.thpt;
            known_results(l).dly = dly';
            known_results(l).coll_p = ra.coll_p';
            known_results(l).succ_p = ra.succ_p';
            known_results(l).n_stas = stas(s);
            known_results(l).agg = aggs(a);
            known_results(l).traffic = trafs(t);
            l=l+1;

        end
    end
end


save("known_results.mat","known_results");
