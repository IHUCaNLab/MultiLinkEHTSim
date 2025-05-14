function [results, sta_table] = single_sim(settings, constants)
%SINGLE_SIM  Runs a SINGLE simulation given the two structures SETTINGS 
%and CONSTANTS.
%
%[results, sta_table] = SINGLE_SIM(settings, constants)
%                       where:
%                       - settings : structure containing the settings.
%                       - constants: structure containing the constants
%                       (see consts.m file).
%                       - results  : simulation results.
%                       - sta_table: final state of the stations' table.

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


% just to make typing more lazy...
s = settings;   c = constants;

n_stas = length(s.sta_exes);

% preallocating table of stas
sta_table = zeros(n_stas,17);

shadowing = ones(n_stas, 1) .* -1;
shadowing_times = ones(n_stas, 1) .* -inf;
nwalls = zeros(n_stas,1);

% initialize every station
for i=1:n_stas
    sta_table(i,1) = i;                     % station id
    sta_table(i,2) = randi([0 c.cw_min]);   % bo counter
    sta_table(i,3) = c.cw_min;              % current contention window size
    
    sta_table(i,4) = 0;                 % number of tx successes
    sta_table(i,5) = 0;                 % number of tx collisions
    sta_table(i,6) = 0;                 % transmitted data
    sta_table(i,7) = 0;                 % tx throughput
    sta_table(i,8) = 0;                 % tx tries
    sta_table(i,9) = 0;                 % dropped packets
    sta_table(i,10) = 0;                % retransmission number

    % STAs positions
    sta_table(i, 11) = s.sta_exes(i);       % x position
    sta_table(i, 12) = s.sta_yis(i);        % y position
    sta_table(i, 13) = 1.2;                 % fixed height
    

	dice = rand();
	
	if dice <= 0.2
		nwalls(i) = 0;
	elseif dice > 0.2 && dice <= 0.6
		nwalls(i) = 1;
	else
		nwalls(i) = 2;
	end

	shadowing(i) = 5 * randn();

    sta_table(i,14) = -1;       % time of last channel realization

    sta_table(i,15) = 0;        % delivered packets on 2.4
    sta_table(i,16) = 0;        % delivered packets on 5
    sta_table(i,17) = 0;        % delivered packets on 6

    tx_dly{i} = [];             % structure for tx delays
    acc_dly{i} = [];            % structure for access delays
    used_mcs{i} = [];           % structure for used MCSs
    tx_times_store{i} = [];

    dataframes_plus_control{i}= [];
    dataframes_times{i}= [];
    dataframes_plus_control_no_ifs{i}= [];
end

%%%%%%%%%%%%%%%%%

for h=1:14
    mcs_vs_bits{h} = nan(0,3);
end



dropped_rtx_limit = zeros(1,n_stas);
dropped_lifetime = zeros(1,n_stas);


% only for debug purposes... plots the scenario
if s.show_scenario
    figure (1);
    plot(s.ap_x, s.ap_y, 'r*','LineWidth',2);
    hold on
    plot(sta_table(:,11), sta_table(:,12),'bx','LineWidth',1);
end


% creating a structure where to store the channel responses
for j=1:3
    for k=1:n_stas
        ch_responses{j}.node{k} = ones(bw2nsd(s.ch_bws(j)),c.spatial_streams*c.spatial_streams).*-1;
        prev_mcs{k} = [-1 -1 -1];
    end
end



% creating tx buffers for all stas
% poisson arrivals...

% converting from traffic rate to arrival rate
if s.max_pkt_size == s.min_pkt_size
    lambdas = s.traffic_rates./(s.max_pkt_size);
    equal_size = 1;
else
    avg_pkt_size = round(mean(s.max_pkt_size,s.min_pkt_size));
    lambdas = s.traffic_rates./(avg_pkt_size);
    equal_size = 0;
end

t_interval = 0.05;

for k=1:n_stas
    pkts = [];
    npktgen = round(lambdas(k) * s.sim_time * 1.15);
    if equal_size == 1
        iarrtimes = exprnd(1/lambdas(k), npktgen, 1);  
        pkts = ones(npktgen,1)*s.min_pkt_size;
        arrival_times = cumsum(iarrtimes,1);
    else
        pkts = randi([s.min_pkt_size, s.max_pkt_size],npktgen,1);
        tempt = (pkts)/s.traffic_rates(k);
        iarrtimes = tempt .* (-log(rand(npktgen,1)));
        arrival_times = cumsum(iarrtimes,1);
        while arrival_times(end) < (s.sim_time)
            pkts = [pkts; randi([s.min_pkt_size, s.max_pkt_size],500,1)];
            tempt = (pkts)/s.traffic_rates(k);
            iarrtimes = tempt .* (-log(rand(length(tempt),1)));
            arrival_times = cumsum(iarrtimes,1);
        end

    end

    

    mask = arrival_times(:) <= s.sim_time;
    traffic{k}(:,1) = 1:1:nnz(mask);
    traffic{k}(:,2) = arrival_times(mask);
    traffic{k}(:,3:4) = zeros(nnz(mask),2);
    % lifetime
    traffic{k}(:,5) = traffic{k}(:,2) + s.mpdu_lifetime;

    traffic{k}(:,6) = ones(nnz(mask),1) * -1;       % time reaches hoq

    traffic{k}(:,7) = pkts(mask);

    traffic{k} = single(traffic{k});
end

for k=1:n_stas
    pointer(k) = binsearch(traffic{k}(:,2), t_interval*2);
    intbuf{k} = traffic{k}(1:pointer(k),:);
end

tx_traces = traffic;

% preallocating space for the RX traces
for ll=1:n_stas
    rx_traces{ll} = zeros(size(traffic{ll},1), size(traffic{ll},2)+1);
    dropped{ll} = [];
end

%%%%%%%% PREALLOCATIONS FOR SPEED %%%%%%%%%%%%%%%%

cur_buf_status = nan(n_stas,1);
pot_tx_times = nan(n_stas,1);
starting_times = nan(n_stas,1);
mindexes = ones(1,length(c.freqs))*-1;
snrs_eff = ones(1,length(c.freqs))*-1;



c_time = 0;     % time counter

rx_idxs = ones(1, n_stas);       % index array for the RX buffers

ch_acc_starting_t = ones(1,n_stas) * -1;

t_th = t_interval;

while(c_time < s.sim_time)

    % getting the current buffer status for all stations
    for k=1:n_stas
        %cur_buf_status(k) = nnz(traffic{k}(:,2) <= c_time);
        if ~isempty(intbuf{k})
            cur_buf_status(k) =binsearch(intbuf{k}(:,2), c_time);
        else
            cur_buf_status(k)=0;
        end
    end

    % computing the potential transmission times
    for k=1:n_stas
        %if cur_buf_status(k) > 0
        if cur_buf_status(k) >= 3
            pot_tx_times(k) = c_time + sta_table(k,2) * c.slot_size;
            starting_times(k) = c_time;

            if ch_acc_starting_t(k) < 0
                ch_acc_starting_t(k) = c_time;
            end
        else
            if size(intbuf{k},1) >= 3
                pot_tx_times(k) = intbuf{k}(3,2) + sta_table(k,2) * c.slot_size;
                starting_times(k) = intbuf{k}(3,2);
    
                if ch_acc_starting_t(k) < 0
                    ch_acc_starting_t(k) = intbuf{k}(3,2);
                end
            else
                pot_tx_times(k) = inf + sta_table(k,2) * c.slot_size;
                starting_times(k) = inf;
    
                if ch_acc_starting_t(k) < 0
                    ch_acc_starting_t(k) = inf;
                end
            end

        end
        
    end

    if all(pot_tx_times == inf)
        break;
    end

    % getting the stations that try to transmit in the same time +/- 2 micros (2 micros is the
    % average time for 802.11 chipsets to switch from RX to TX)
    mask = pot_tx_times >= (min(pot_tx_times) - 2e-6) & pot_tx_times <= (min(pot_tx_times) + 2e-6);

    % if only 1 sta tries in that time frame, no collision
    if nnz(mask) == 1
        % no collision
        % finding the winner
        [~, winner] = min(pot_tx_times);

        % incrementing the number of transmission tries
        sta_table(winner,8) = sta_table(winner,8) +1;

        % incrementing the number of channel access successes
        sta_table(winner,4) = sta_table(winner,4) +1;

        % resetting the contention window and collision counter
        sta_table(winner,3) = c.cw_min;
        sta_table(winner, 10) = 0;
        
        % computing how much should be subtracted from all other backoff
        % counters...
        % NOTE: the rounding is needed because in some cases, due to not
        % perfect sync, a small error (in the e-21 order of magnitude) can
        % appear and propagate due to the division by the slot size and the
        % consequent flooring.
        to_sub = floor(round((min(pot_tx_times)-starting_times),10)/c.slot_size);
        % ...but stations that didn't have any traffic before a lot of
        % time, we don't subtract anything.
        to_sub(to_sub < 0) = 0;
        sta_table(:,2) = sta_table(:,2) - to_sub;

        % new value for backoff counter for the winner
        sta_table(winner,2) = randi([0 sta_table(winner, 3)]);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HERE THE RADIO PART BEGINS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rates = ones(1,length(c.freqs))*-1;
        
        if (c_time - shadowing_times(winner)) > s.shad_refresh
            shadowing(winner) = 5*randn();        % 5 dbs of variance
            shadowing_times(winner) = c_time;
        end
			
		if (c_time - sta_table(winner,14) > s.ch_refresh) || (c_time - shadowing_times(winner) > s.shad_refresh)
		
			for i=1:length(c.freqs)
                nsd = bw2nsd(s.ch_bws(i));
				[rho, rx_p] = get_rho(c.sta_tx_p, sta_table(winner,11),s.ap_x, sta_table(winner,12), s.ap_y, sta_table(winner,13), s.ap_z, c.sta_ant_gain, c.ap_ant_gain, c.freqs(i), shadowing(winner), nwalls(winner), c);
				rho_db(i) = 10 * log10(rho);
				rx_p_db(i) = rx_p;
			end
			
			while (any(rx_p_db < s.cca_below_lims))
				shadowing(winner) = 5*randn();        % 5 dbs of variance
				shadowing_times(winner) = c_time;
				for i=1:length(c.freqs)
                    nsd = bw2nsd(s.ch_bws(i));
					[rho, rx_p] = get_rho(c.sta_tx_p, sta_table(winner,11),s.ap_x, sta_table(winner,12), s.ap_y, sta_table(winner,13), s.ap_z, c.sta_ant_gain, c.ap_ant_gain, c.freqs(i), shadowing(winner), nwalls(winner), c);
					rho_db(i) = 10 * log10(rho);
					rx_p_db(i) = rx_p;
				end
			end
			
            for i=1:length(c.freqs)
                nsd = bw2nsd(s.ch_bws(i));
                if s.model == "Mod_C"
                    snrs = get_snr_mod_c(rho_db(i), nsd, c.spatial_streams, c);
                elseif s.model == "Stat_mod"
                    snrs = get_snr_cdf(rho_db(i), nsd, c.spatial_streams, c);
                else
                    error("Error: invalid channel model.\n");
                end
                    
			    mcs_index = get_mcs(rx_p_db(i), s.ch_bws(i));
			    snrs_eff(i) = compute_snreff(snrs, mcs_index, c);
			    rr(i) = snrs_eff(i);
			    mindexes(i) = mcs_index;
			    prev_mcs{winner}(i) = mcs_index;
    
			    rates(i) = mcs2rate(mcs_index, c.spatial_streams, s.ch_bws(i), c.guard_int, c);
			    mcsar(i) = mcs_index;
            end
			
        else
            for i=1:length(c.freqs)
			    rates(i) = mcs2rate(prev_mcs{winner}(i), c.spatial_streams, s.ch_bws(i), c.guard_int, c);
			    mcsar(i) = prev_mcs{winner}(i);
            end
		end


        if (c_time - sta_table(winner, 14)) > s.ch_refresh
            sta_table(winner, 14) = c_time;
        end

        rates = rates * 1e6;

        % storing the used MCS indexes
        used_mcs{winner} = [used_mcs{winner}; mindexes];

        n_pkts_in_buffer =binsearch(intbuf{winner}(:,2), pot_tx_times(winner));
        n_pkts = min(n_pkts_in_buffer, s.num_agg_pkts);

        packets = intbuf{winner}(1:n_pkts,7);

        [PP,tx_times] = aggregate(s.agg_meth, packets, rates,3.008e-3,s.ch_bws,c);

        Ts = c.rts + c.cts + (3 * c.sifs) + max(tx_times) + c.ack + c.difs;
        

        for p=1:length(packets)     
            if intbuf{winner}(p, 6) == -1
                intbuf{winner}(p, 6) = pot_tx_times(winner);
            end
        end

        if s.perfect_rx == 0
            received_mask=[];
            for i=1:length(c.freqs)

                if ~isempty(PP{i}) 
                    threshold = snreff2per(snrs_eff(i), mindexes(i), c);
                    random_eval = rand(length(PP{i}), 1);
                    received{i} = random_eval >= threshold;
                    received_mask = [received_mask; random_eval >= threshold];
                else
                    received{i} = 0;
                end
            end
            received_mask = logical(received_mask);
            % if none of the packets were correctly delivered (VERY RARE, BUT
            % YOU NEVER KNOW AND BETTER COVER ALL POSSIBILITIES)
            if nnz(received_mask) == 0
                % not sending any ack at all
                Ts = c.rts + c.cts + (3 * c.sifs) + max(tx_times) + c.difs;
            end
        else
            received_mask = ones(length(packets),1);        % all packets correctly received
            received_mask = logical(received_mask);

            for i=1:length(c.freqs)
                received{i} = ones(length(PP{i}),1);        
            end
        end

        for l=1:14
            ssizes(l) = size(mcs_vs_bits{l},1);
        end

        for l=1:3
            mcs_vs_bits{mcsar(l)+1}(ssizes(mcsar(l)+1)+1,l) = sum(packets(PP{l}) );
        end


        times_to_store = tx_times;
        times_to_store(1) = Ts;        

        dataframes_plus_control{winner} = [dataframes_plus_control{winner}; times_to_store];
        times_to_store(1) = times_to_store(1) - (3* c.sifs + c.difs);
        dataframes_plus_control_no_ifs{winner} = [dataframes_plus_control_no_ifs{winner}; times_to_store];
        dataframes_times{winner} = [dataframes_times{winner}; tx_times];

        c_time = pot_tx_times(winner) + Ts;

        % incrementing how many packets have been delivered for each band
        sta_table(winner, 15) = sta_table(winner, 15) + nnz(received{1});      % 2.4 GHz
        sta_table(winner, 16) = sta_table(winner, 16) + nnz(received{2});      % 5 GHz
        sta_table(winner, 17) = sta_table(winner, 17) + nnz(received{3});      % 6 GHz

        tried_pkts = intbuf{winner}(1:length(packets),:);
        for j=1:3
            tried_pkts(PP{j}, 8) = j;
        end
        dlvrd_pkts = tried_pkts(received_mask,:);       % select the packets that were delivered
        dlvrd_pkts(:,4) = c_time;        % store the arrival time
        
        
        rx_traces{winner}(rx_idxs(winner):rx_idxs(winner)+nnz(received_mask)-1,:) = dlvrd_pkts;     % store them to the corresponding RX trace...
        
        rx_idxs(winner) = rx_idxs(winner) + nnz(received_mask);     % increment the corresponding shift

        mask_to_inc = not(received_mask);
        intbuf{winner}(mask_to_inc, 3) = intbuf{winner}(mask_to_inc, 3) +1;     % incrementing retry number for failed pkts
        mask_to_del = received_mask;

        intbuf{winner}(mask_to_del, :) = [];      % remove delivered pkts from TX buffer

        mask_del_d = intbuf{winner}(:,3) > s.pkt_max_retx;

        if nnz(mask_del_d) > 0
            if isempty(dropped{winner})
                dropped{winner}(1:nnz(mask_del_d),:) = intbuf{winner}(mask_del_d,:);
            else
                dropped{winner}(end+1:end+nnz(mask_del_d),:) = intbuf{winner}(mask_del_d,:);
            end
        end

        dropped_rtx_limit(winner) = dropped_rtx_limit(winner) + nnz(mask_del_d); 
        intbuf{winner}(intbuf{winner}(:,3) > s.pkt_max_retx , :) = [];      % removing pkts that were transmitted too many times

        if s.mpdu_lifetime ~= inf
            for k=1:n_stas
                mask_expired = c_time > intbuf{k}(:,5);
                intbuf{k}(mask_expired,:) = [];
                dropped_lifetime(k) = dropped_lifetime(k) + nnz(mask_expired);
            end
        end
  
        tx_dly{winner} = [tx_dly{winner}; (pot_tx_times(winner) - ch_acc_starting_t(winner)) + Ts];
        acc_dly{winner} = [acc_dly{winner}; (pot_tx_times(winner) - ch_acc_starting_t(winner))];
        ch_acc_starting_t(winner) = -1;

    else
        % collision
        colliding = sta_table(mask,1);

        others = not(mask);

        % incrementing number of tries for the colliding nodes
        sta_table(colliding,8) = sta_table(colliding,8) +1;
        % increment number of station collisions for colliding nodes
        sta_table(colliding,5) = sta_table(colliding,5) +1;

        % for each colliding station...
        for sidx = colliding'
            % determine the number of packets currently in buffer 
            n_pkts_in_buffer =binsearch(intbuf{sidx}(:,2), pot_tx_times(sidx));
            % if there are enough packets in the buffer, transmit the maximum
            % number of packets aggregable, otherwise transmit whatever we have
            % in buffer.
            n_pkts = min(n_pkts_in_buffer, s.num_agg_pkts);


            for p=1:n_pkts
                if intbuf{sidx}(p, 6) == -1
                    intbuf{sidx}(p, 6) = pot_tx_times(sidx);
                end
            end

            % mask to select the packets that still have some retrial at their
            % disposal...
            mask_retry = intbuf{sidx}(1:n_pkts,3) < c.max_retrial;
            % and increment their retrial counter.
            intbuf{sidx}(mask_retry , 3) = intbuf{sidx}(mask_retry , 3) +1;

            % mask to select the other packets (the ones that exceeded the retry limit)
            mask_dropped = not(mask_retry);

            % remove the expired packets
            if isempty(dropped{sidx})
                if nnz(mask_dropped)>0
                    dropped{sidx}(1:nnz(mask_dropped),:) = intbuf{sidx}(mask_dropped,:);
                end
            else
                if nnz(mask_dropped)>0
                    dropped{sidx}(end+1:end+nnz(mask_dropped),:) = intbuf{sidx}(mask_dropped,:);
                end
            end

            intbuf{sidx}(mask_dropped, :) = [];

            % increment the number of dropped packets
            sta_table(sidx, 9) = sta_table(sidx, 9) + nnz(mask_dropped);

            dropped_rtx_limit(sidx) = dropped_rtx_limit(sidx) + nnz(mask_dropped);

            cur_coll_cntr = sta_table(sidx, 10);
            % if the station global collision counter is less than the number of max retries,
            % increment it...
            if cur_coll_cntr < c.max_retrial
                sta_table(sidx, 3) = min((2^(cur_coll_cntr+1)*(c.cw_min+1)-1), c.cw_max);
                sta_table(sidx, 10) = sta_table(sidx, 10) + 1;      % incrementing number of contention window increments
            else
                sta_table(sidx, 3) = c.cw_min;                      % resetting the contention window to cw min
                sta_table(sidx, 10) = 0;                            % resetting number of contention window increments
            end
        end

        % computing how much should be subtracted from all other backoff
        % counters...
        % NOTE: the rounding is needed because in some cases, due to not
        % perfect sync, a small error (in the e-21 order of magnitude) can
        % appear and propagate due to the division by the slot size and the
        % consequent flooring.
        to_sub = floor(round((min(pot_tx_times)-starting_times),10)/c.slot_size);
        
        mask_not_competing = to_sub < 0;
        % ...but stations that didn't have any traffic before a lot of
        % time, we don't subtract anything. (-1 gets cancelled later on)
        to_sub(mask_not_competing) = -1;
    
        sta_table(:,2) = sta_table(:,2) - (to_sub+1);

        mask_others_competing = others & not(mask_not_competing);
        sta_table(mask_others_competing, 2) = sta_table(mask_others_competing, 2) -1;

        % new backoff counters for the colliding stas
        for cidx = colliding'
            sta_table(cidx,2) = randi([0 sta_table(cidx,3)]);
        end

        ch_acc_starting_t(colliding) = -1;
        % incrementing time
        c_time = pot_tx_times(colliding(1)) + c.Tc;

        if s.mpdu_lifetime ~= inf
            for k=1:n_stas
                mask_expired = c_time > intbuf{k}(:,5);
                dropped_lifetime(k) = dropped_lifetime(k) + nnz(mask_expired);
                intbuf{k}(mask_expired,:) = [];
            end
        end
    end

    if (c_time >= t_th)
        for k=1:n_stas
            ntoadd = min(round(lambdas(k) * 2*t_interval), size(traffic{k},1)-pointer(k));
            intbuf{k} = [intbuf{k}; traffic{k}(pointer(k)+1 : pointer(k)+ntoadd,:)];
            pointer(k) = pointer(k) + ntoadd;
        end
        t_th = t_th + t_interval;
    end
end

num_succ = sta_table(:,4);
num_tries = sta_table(:,8);
num_colls = sta_table(:,5);


%%%%%%%%%%%    STATISTICS        %%%%%%%%%%%%%%%%%%

% computing throughput per station
for k=1:n_stas
    results.thpt(k) = sum(rx_traces{k}(:,7)) / s.sim_time;
    agg_thpt_24(k) = sum(rx_traces{k}(rx_traces{k}(:,8) == 1, 7))/s.sim_time;
    agg_thpt_5(k) = sum(rx_traces{k}(rx_traces{k}(:,8) == 2, 7))/s.sim_time;
    agg_thpt_6(k) = sum(rx_traces{k}(rx_traces{k}(:,8) == 3, 7))/s.sim_time;
end

% computing aggregated throughput per band

results.agg_thpt_24 = sum(agg_thpt_24);
results.agg_thpt_5 = sum(agg_thpt_5);
results.agg_thpt_6 = sum(agg_thpt_6);


% computing mean throughput per band
results.mean_thpt_24 = results.agg_thpt_24 / n_stas;
results.mean_thpt_5 = results.agg_thpt_5 / n_stas;
results.mean_thpt_6 = results.agg_thpt_6 / n_stas;

results.mcs_vs_bits = mcs_vs_bits;

results.succ_p = num_succ ./ num_tries;     % success probabilities
results.coll_p = num_colls ./ num_tries;    % collision probabilities
results.access_dly = acc_dly;
results.tx_dly = tx_dly;
results.agg_thpt = sum(results.thpt);

results.fairness = (sum(results.thpt)^2) / (n_stas * sum(results.thpt .^ 2));

results.dropped_lifetime = dropped_lifetime;
results.dropped_rtx_limit = dropped_rtx_limit;
results.succ_pkts = sum(sta_table(:,15:17),2);

results.sta_table = sta_table;

% optimizing full results storage...

for k=1:n_stas
    results.dly_perc_95(k) = single(prctile(rx_traces{k}(:,4) - rx_traces{k}(:,2), 95));
    results.dly_perc_99(k) = single(prctile(rx_traces{k}(:,4) - rx_traces{k}(:,2), 99));

    dly{k} = single(rx_traces{k}(1:rx_idxs(k)-1,4) - rx_traces{k}(1:rx_idxs(k)-1,2));

    hoq_dly{k} = single(rx_traces{k}(1:rx_idxs(k)-1,4) - rx_traces{k}(1:rx_idxs(k)-1,6));

    ch_occ_dataframes_per_sta(k,:) = sum(dataframes_times{k});
    avg_df_times(k,:) = mean(dataframes_times{k},1);

    ch_occ_df_plus_control_per_sta(k,:) = sum(dataframes_plus_control{k});
    avg_df_times_control(k,:) = mean(dataframes_plus_control{k},1);

    ch_occ_df_plus_control_no_ifs_per_sta(k,:) = sum(dataframes_plus_control_no_ifs{k});
    avg_df_times_control_no_ifs(k,:) = mean(dataframes_plus_control_no_ifs{k},1);
end

results.ch_occ_dataframes = ch_occ_dataframes_per_sta ./ s.sim_time;
results.ch_occ_df_plus_control = ch_occ_df_plus_control_per_sta ./ s.sim_time;
results.ch_occ_df_plus_control_no_ifs = ch_occ_df_plus_control_no_ifs_per_sta ./ s.sim_time;

results.dataframes_times = dataframes_plus_control;
results.df_plus_control = dataframes_plus_control;
results.df_plus_control_no_ifs = dataframes_plus_control_no_ifs;
results.nwalls = nwalls;

results.avgdataframes_times = avg_df_times;
results.avgdf_plus_control = avg_df_times_control;
results.avgdf_plus_control_no_ifs = avg_df_times_control_no_ifs;

results.dly_per_sta = dly;
results.hoq_dly_per_sta = hoq_dly;

results.dly_perc_95_mean = mean(results.dly_perc_95);
results.dly_perc_95_mean = mean(results.dly_perc_95);

end
