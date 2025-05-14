%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrity_test.m
% 
% Performs an integrity test for the main script.
% A simulation is performed with a fixed seed for the random generator and
% the results are compared with known good results.
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

% known results
load known_results.mat


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
settings.agg_meth = 1;
settings.ap_x = 10;
settings.ap_y = 10;
settings.min_pkt_size = 1500*8;
settings.max_pkt_size = 1500*8;
settings.use_traces = 0;
settings.model = "Stat_mod";

n_passed=0;

for i=1:length(known_results)
    rng(4707);

    nstas = known_results(i).n_stas;
    settings.traffic_rates= ones(1,nstas) * known_results(i).traffic;
    settings.num_agg_pkts = known_results(i).agg;
    rand_rs = 8 * sqrt(rand(nstas,1));
    rand_ths = 2*pi*rand(nstas,1);
    [xes,yis] = pol2cart(rand_ths, rand_rs);
    settings.sta_exes = settings.ap_x + xes;
    settings.sta_yis = settings.ap_y + yis;

    [ra, stat] = single_sim(settings,c);

    dly = [];
    for k=1:nstas
        dly=[dly; mean(ra.dly_per_sta{k})];
    end
    dly=round(dly,4);

    if max(abs(ra.thpt - known_results(i).thpt)) == 0
        if max(abs(round(dly' - known_results(i).dly,6)))==0
            if max(abs(round(ra.coll_p' - known_results(i).coll_p,6)))==0
                if max(abs(round(ra.succ_p' - known_results(i).succ_p,6)))==0
                    fprintf("TEST CASE %d PASSED\n",i);
                    n_passed = n_passed +1;
                else
                    fprintf("TEST CASE %d FAILED\n",i);
                    fprintf(" - Succ. prob.\n")
                    break;
                end
            else
                fprintf("TEST CASE %d FAILED\n",i);
                fprintf(" - Coll. prob.\n")
                break;
            end
        else
            fprintf("TEST CASE %d FAILED\n",i);
            fprintf(" - Dly\n")
            break;
        end
    else
        fprintf("TEST CASE %d FAILED\n",i);
        fprintf(" - Thpt\n")
        break;
    end
end

if n_passed == length(known_results)
    fprintf("=============================\n");
    fprintf("ALL TESTS PASSED.\n");
else
    fprintf("=============================\n");
    fprintf("TESTING FAILED.\n");
end
