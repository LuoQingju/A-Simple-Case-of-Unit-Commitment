%% A Simple Case of Unit Commitment (UC)

% yalmip modeling, gurobi solving

% Dependency libraries: yalmip, gurobi, matpower

% Author: Qingju Luo
% Email: luoqingju@qq.com

% The author's ability is limited, there will inevitably be errors and inadequacies, please criticize and correct!

% For learning and communication only!!!
% For learning and communication only!!!
% For learning and communication only!!!

clc
clear
define_constants;

mpc = case9;

%% Get Data

baseMVA = mpc.baseMVA;
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;

nb = size(bus, 1); %% number of buses
nl = size(branch, 1); %% number of branches
ng = size(gen, 1); %% number of dispatchable injections

np = 24; %% number of periods

load_rate = 1 + 0.95 * cos(0:4*pi/(np - 1):4*pi); % load rate

gen(:, RAMP_AGC) = gen(:, PMAX) * 0.01; % ramp rate for load following/AGC (MW/min)
ramp_rate = 60 .* gen(:, RAMP_AGC) / baseMVA; % ramp rate for load following/AGC (p.u./h)

RTmin = [3; 3; 4]; % Minimum runtime
DTmin = [7; 6; 8]; % Minimum downtime

Sg0 = gen(:, GEN_STATUS); % status, > 0 - in service, <= 0 - out of service
Sg0(Sg0 > 0) = 1;
Sg0(Sg0 <= 0) = 0;

Pg0 = gen(:, PG) / baseMVA; % real power output (p.u.)
Pg0 = Pg0 .* Sg0;


Pmax = gen(:, PMAX) / baseMVA; % maximum real power output (p.u.)
Pmin = gen(:, PMIN) / baseMVA; % minimum real power output (p.u.)

% generator costs
Cup = gencost(:, STARTUP); % startup
Cdown = gencost(:, SHUTDOWN); % shutdown
C2 = gencost(:, COST); % quadratic
C1 = gencost(:, COST+1); % linear
C0 = gencost(:, COST+2); % constant

Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %% connection matrix for generators & buses


ref = find(bus(:, BUS_TYPE) == REF); % reference bus
Va_ref = bus(ref, VA) * (pi / 180); % reference bus voltage angle (radians)

Pd = bus(:, PD) / baseMVA; % real power demand (p.u.)
Gs = bus(:, GS) / baseMVA; % shunt conductance (p.u. demanded at V = 1.0 p.u.)

% Builds the B matrices and phase shift injections for DC power flow.
[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);

flow_max = branch(:, RATE_A) / baseMVA; % branch flow limit (p.u.)
il = find(flow_max ~= 0); % flow constrained lines
if ~isempty(il)
    upf = flow_max(il) - Pfinj(il);
    upt = flow_max(il) + Pfinj(il);
end

%% Define variables
Sg = binvar(ng, np, 'full'); % generator state
Up = binvar(ng, np, 'full'); % startup
Dn = binvar(ng, np, 'full'); % shutdown
Pg = sdpvar(ng, np, 'full'); % real power output (p.u.)
Va = sdpvar(nb, np, 'full'); % voltage angle (radians)

%% Define the constraints and objective function
con = [];
obj = 0;

for t = 1:np

    %% real power balance
    con = con + [Cg * Pg(:, t) == B * Va(:, t) + Pbusinj + Gs + Pd .* load_rate(t)]; %#ok<*NBRAK>

    %% reference bus voltage angle constraint
    con = con + [Va(ref, t) == Va_ref];

    %% flow limit
    if ~isempty(il)
        con = con + [-upt <= Bf(il, :) * Va(:, t) <= upf]; %#ok<*CHAIN>
    end

    %% power output limit
    con = con + [Pmin .* Sg(:, t) <= Pg(:, t)];
    con = con + [Pg(:, t) <= Pmax .* Sg(:, t)];

    %% ramp limit
    if t > 1
        con = con + [-ramp_rate <= Pg(:, t) - Pg(:, t-1) <= ramp_rate];
    else
        con = con + [-ramp_rate <= Pg(:, t) - Pg0 <= ramp_rate];
    end

    %% logical constraint
    if t > 1
        con = con + [Up(:, t) - Dn(:, t) == Sg(:, t) - Sg(:, t-1)];
    else
        con = con + [Up(:, t) - Dn(:, t) == Sg(:, t) - Sg0];
    end

    %% minimum runtime and downtime constraints
    for i = 1:ng
        RTidx = t : -1 : (t - RTmin(i) + 1);
        DTidx = t : -1 : (t - DTmin(i) + 1);
        RTidx(RTidx < 1) = [];
        DTidx(DTidx < 1) = [];

        con = con + [sum(Up(i, RTidx), 2) <= Sg(i, t)];
        con = con + [sum(Dn(i, DTidx), 2) <= (1 - Sg(i, t))];
    end

    %% objective function
    obj = obj + Pg(:, t)' * (C2 .* Pg(:, t)) + C1' * Pg(:, t) + C0' * Sg(:, t) + ...
        Cup' * Up(:, t) + Cdown' * Dn(:, t);

end

%% Solve
ops = sdpsettings('solver', 'gurobi', 'verbose', 2, 'savesolveroutput', 1);
sol = solvesdp(con, obj, ops);

if sol.problem ~= 0
    disp(sol.info)
    return
end

obj = value(obj);

Sg = value(Sg);
Up = value(Up);
Dn = value(Dn);
Pg = value(Pg);
Va = value(Va);
