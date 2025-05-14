function [PP, t] = aggregateTBA(P,r,rho,c)
%AGGREGATETBA  Aggregates packets so to balance the time duration of
%the three TXOPs.
%
% [PP, t] = AGGREGATETBA(P,r,tau,rho)
%         where:
%         - P: vector of packet sizes [bits]
%         - r: vector of rates sorted [bps]
%         - rho: TXOP limit [s]

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

[Psrtd, pidx] = sort(P,'descend');
[rsrtd, ridx] = sort(r,'descend');

Amax = 1048575 *8;   % maximum number of bits that fit an A-MPDU
nu = length(Psrtd);

tPP{1} = nan(nu,1);
tPP{2} = nan(nu,1);
tPP{3} = nan(nu,1);
PP{1} = [];
PP{2} = [];
PP{3} = [];
PPu = nan(nu,1);
tPP_idx = [1,1,1];

L = [c.PLCP_head, c.PLCP_head, c.PLCP_head];
t = c.T_preamble + c.PLCP_head./rsrtd;

W = [(Psrtd + c.MAC_head)/rsrtd(1), (Psrtd + c.MAC_head)/rsrtd(2), (Psrtd + c.MAC_head)/rsrtd(3)];

tau = ((rsrtd(3)*((nu * c.MAC_head + c.PLCP_head) + sum(Psrtd)))/(sum(rsrtd))) / rsrtd(3);

ppu_idx=1;
for i=1:nu
    
    unass = 1;   % flag to signal if a packet is left unassigned
    for j=1:3
        if (t(j)+W(i,j)<= tau) && (L(j) + c.MAC_head + Psrtd(i) <= Amax)
            t(j) = t(j) + W(i,j);
            tPP{j}(tPP_idx(j)) = i;
            tPP_idx(j) = tPP_idx(j) +1;
            L(j) = L(j) + c.MAC_head + Psrtd(i);
            unass = 0;
            break;
        end
    end

    if unass == 1
        PPu(ppu_idx) = i;
        ppu_idx = ppu_idx +1;
    end
end

for u=1:ppu_idx-1
    [~, minind] = min(t);
    if (t(minind) + W(PPu(u),minind) <= rho) && (L(minind)+c.MAC_head + Psrtd(PPu(u))<=Amax)
        tPP{minind}(tPP_idx(minind)) = PPu(u);
        tPP_idx(minind) = tPP_idx(minind)+1;
        t(minind) = t(minind) + W(PPu(u),minind);
        L(minind) = L(minind) + c.MAC_head + Psrtd(PPu(u));
    end
end

for j=1:3
    PP{ridx(j)} = pidx(tPP{j}(1:tPP_idx(j)-1));
end

end