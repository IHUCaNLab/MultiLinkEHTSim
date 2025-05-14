function [PP,t] = aggregate(method, P, r, rho, bws,c)
%AGGREGATE  Aggregates packets according to a specific method.
%
%[PP, t] = AGGREGATE(method, P, r, rho, bws, c)
%          where:
%          - method: aggregation method. 1) TBA; 2) ESA; 3)BPA.
%          - P     : MPDUs to be aggregated.
%          - r     : rates vector
%          - rho   : TXOP limit
%          - bws   : channel bandwidths vector
%          - c     : common constants structure
%          - PP    : structure containing the agggregated MPDU indexes
%          - t     : vector containing the tx durations for the three
%          channels

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





N = length(P);
t = [NaN,NaN,NaN];

if N == 1
    [~, fi] = max(r);
    PP{1} = [];
    PP{2} = [];
    PP{3} = [];
    PP{fi} = [PP{fi}; 1];
elseif N == 2
    [~,fi] = max(r);
    [~,si] = min(r);
    bb = 1:3;
    mid = bb(bb ~= fi & bb ~= si);
    PP{fi} = [PP{fi}; 1];
    PP{mid} = [PP{mid}; 2];
    PP{si} = [];
else
    switch method
        case 1
            % Time Balanced Aggregation
            [PP, t] = aggregateTBA(P,r,rho,c);

        case 2
            % Equal Split Aggregation
            NperLink = round(N/3);
            PP{1} = 1:1:NperLink;
            PP{2} = NperLink+1:1:NperLink*2;
            PP{3} = (NperLink*2)+1:1:N;
 
        case 3
            % Bandwidth Proportional Aggregation
            beta = bws ./ min(bws);
            jk = N/sum(beta);
            Ml = round(jk * beta);
            PP{1} = 1:1:Ml(1);
            PP{2} = Ml(1)+1 : 1 : Ml(1)+Ml(2);
            PP{3} = Ml(1)+Ml(2)+1:1:N;

        otherwise
            error("Aggregation method non existent!");
    end
end

    % TBA already returns the vector t, so we need to compute it only for
    % the other methods. This does the trick.
    if all(isnan(t))

        for j=1:3
            t(j) = c.T_preamble + c.PLCP_head/r(j) + (length(PP{j})*c.MAC_head + sum(P(PP{j})))/r(j);
        end
    end

end
