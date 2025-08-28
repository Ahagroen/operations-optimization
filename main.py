import gurobipy as gp
from gurobipy import GRB, quicksum

NF = ["F1", "F2", "F3", "F4"]           
K  = ["AC1", "AC2"]           # aircraft
MT = ["M1", "M2"]           # maintenance stations
A  = ["A1", "A2", "A3"]           # airports
o, t = "o", "t"      # dummy source/sink
KT   = 2                                    # total number of aircraft used (if needed for V computing)
TRT  = 45                                    # turn-around time
Tmax = 40*60                                    # max flying time between successive maint ops
DT   = {"F1": 0, "F2": 480, "F3": 600, "F4":720}                   # departure time of flight leg i
AT   = {"F1": 225, "F2": 600, "F3": 710, "F4":825}                   # arrival time of flight leg i
Vmax = int(round(sum([DT[i] for i in NF]) / (Tmax * KT)))           # int
Vset = list(range(1, Vmax+1))
NF_i = NF + [o]      # i-domain where i can be o
NF_j = NF + [t]      # j-domain where j can be t
Oi_a = {(i,a): ... for i in NF for a in A}    # O_{i a} (origin indicator)
Di_a = {(i,a): ... for i in NF for a in A}    # D_{i a} (destination indicator)
FT   = {i: AT[i]-DT[i] for i in NF}                   # flight duration of leg i (if needed elsewhere)
bij  = {(i,j): ... for i in NF for j in NF}   # through value between legs i->j
Cmax = 15                                    # max number of take-offs between successive maint ops
MP   = {m: ... for m in MT}                   # workforce capacity at station m
ET   = {m: ... for m in MT}                   # close time for station m
Mb   = {(m,a): ... for m in MT for a in A}    # 1 if station m located at airport a
MAT  = 8*60                                    # time required to perform maintenance
PC   = {(k,m): 500 for k in K for m in MT}    # penalty cost (as in objective)
Mbig = 10**7                                   # big-M (arbitrarily large per your note)

# Helper: through-value b for arcs involving o/t -> treat as zero unless you define otherwise
def b_ext(i, j):
    if (i in NF) and (j in NF):
        return bij[i, j]
    else:
        return 0.0

# =========================
# Model
# =========================
m = gp.Model("AircraftMaintenanceRouting")

# =========================
# Decision variables
# =========================
# x_{i j k v}: i in NF ∪ {o}, j in NF ∪ {t}
x = m.addVars(NF_i, NF_j, K, Vset, vtype=GRB.BINARY, name="x")

# y_{i m k v}: i in NF ∪ {o}
y = m.addVars(NF_i, MT, K, Vset, vtype=GRB.BINARY, name="y")

# z_{m j k v}: j in NF ∪ {t}
z = m.addVars(MT, NF_j, K, Vset, vtype=GRB.BINARY, name="z")

# RTAM_{k v} >= 0
RTAM = m.addVars(K, Vset, lb=0.0, name="RTAM")

# Positive-part linearization for (y_{imkv} - MP_m)+  (faithful to the written objective)
p = m.addVars(NF_i, MT, K, Vset, lb=0.0, name="pospart")
for i in NF_i:
    for m_ in MT:
        for k in K:
            for v in Vset:
                # p >= y - MP_m
                m.addConstr(p[i,m_,k,v] >= y[i,m_,k,v] - MP[m_],
                            name=f"pp_ge[{i},{m_},{k},{v}]")
                # p >= 0 already via lb
                # Optional tightening: p <= y (keeps it tight when MP[m_] >= 0)
                m.addConstr(p[i,m_,k,v] <= y[i,m_,k,v],
                            name=f"pp_le_y[{i},{m_},{k},{v}]")

# =========================
# Objective (as provided)
# max  sum b_ij x_ijkv  -  sum PC_km * (y_imkv - MP_m)+
# =========================
obj_gain = quicksum(b_ext(i,j) * x[i,j,k,v]
                    for i in NF_i for j in NF_j for k in K for v in Vset)
obj_pen  = quicksum(PC[k,m_] * p[i,m_,k,v]
                    for i in NF_i for m_ in MT for k in K for v in Vset)
m.setObjective(obj_gain - obj_pen, GRB.MAXIMIZE)

# =========================
# Constraints
# =========================

# (2) Each flight leg i ∈ NF is covered exactly once (by x leaving i OR maintenance taken after i)
for i in NF:
    m.addConstr(
        quicksum(x[i,j,k,v] for j in NF for k in K for v in Vset)   # i -> j (j in NF)
      + quicksum(y[i,m_,k,v] for m_ in MT for k in K for v in Vset) # i -> maintenance
      == 1,
      name=f"(2)_cover_once[{i}]"
    )

# (3) Start at origin o for each (k,v)
for k in K:
    for v in Vset:
        m.addConstr(
            quicksum(x[o,j,k,v] for j in NF)              # o -> flight j
          + quicksum(y[o,m_,k,v] for m_ in MT)            # o -> maintenance m
          == 1, name=f"(3)_start[{k},{v}]"
        )

# (4) End at sink t for each (k,v)
for k in K:
    for v in Vset:
        m.addConstr(
            quicksum(x[i,t,k,v] for i in NF)              # i -> t
          + quicksum(z[m_,t,k,v] for m_ in MT)            # m -> t
          == 1, name=f"(4)_end[{k},{v}]"
        )

# (5) Flow conservation at each i ∈ NF for each (k,v)
#    sum_in (from j or from maintenance) = sum_out (to j or to maintenance)
for i in NF:
    for k in K:
        for v in Vset:
            lhs = quicksum(x[j,i,k,v] for j in NF) + quicksum(z[m_,i,k,v] for m_ in MT)
            rhs = quicksum(x[i,j,k,v] for j in NF) + quicksum(y[i,m_,k,v] for m_ in MT)
            m.addConstr(lhs == rhs, name=f"(5)_flow[{i},{k},{v}]")

# (6) For each station m and aircraft k, maintenance entries equal exits (across all v)
for m_ in MT:
    for k in K:
        lhs = quicksum(y[j,m_,k,v] for j in (NF + [o]) for v in Vset)  # include o if allowed to go to m first
        rhs = quicksum(z[m_,j,k,v] for j in (NF) for v in Vset)        # exits to real flights only
        m.addConstr(lhs == rhs, name=f"(6)_maint_balance[{m_},{k}]")

# (7) Time feasibility on flight->flight arcs: AT_i + TRT <= DT_j if x_{i j k v} = 1
for i in NF:
    for j in NF:
        for k in K:
            for v in Vset:
                m.addConstr(
                    AT[i] + TRT - DT[j] <= Mbig * (1 - x[i,j,k,v]),
                    name=f"(7)_time_x[{i},{j},{k},{v}]"
                )

# (8) Flight connectivity feasibility:
#     sum_{k,v} x_{i j k v} <= sum_{a∈A} D_{i a} * O_{j a},  ∀ i,j∈NF
for i in NF:
    for j in NF:
        rhs = sum(Di_a[i,a] * Oi_a[j,a] for a in A)
        m.addConstr(
            quicksum(x[i,j,k,v] for k in K for v in Vset) <= rhs,
            name=f"(8)_conn[{i},{j}]"
        )

# (9) Time feasibility flight->maintenance: AT_i + MAT <= ET_m if y_{i m k v} = 1
for i in NF:
    for m_ in MT:
        for k in K:
            for v in Vset:
                m.addConstr(
                    AT[i] + MAT - ET[m_] <= Mbig * (1 - y[i,m_,k,v]),
                    name=f"(9)_time_y[{i},{m_},{k},{v}]"
                )
# If y[o,m,k,v] is allowed, you can skip time link for i=o (no AT[o]); or set AT[o]=0 if meaningful.

# (10) Maintenance feasibility for y: only if leg i ends where station m is located
#      sum_{k,v} y_{i m k v} <= sum_{a∈A} D_{i a} * Mb_{m a}
for i in NF:
    for m_ in MT:
        rhs = sum(Di_a[i,a] * Mb[m_,a] for a in A)
        m.addConstr(
            quicksum(y[i,m_,k,v] for k in K for v in Vset) <= rhs,
            name=f"(10)_y_loc[{i},{m_}]"
        )

# (11) Maintenance feasibility for z: can depart to j only if station m is at j's origin airport
#      sum_{k,v} z_{m j k v} <= sum_{a∈A} Mb_{m a} * O_{j a}
for m_ in MT:
    for j in NF:
        rhs = sum(Mb[m_,a] * Oi_a[j,a] for a in A)
        m.addConstr(
            quicksum(z[m_,j,k,v] for k in K for v in Vset) <= rhs,
            name=f"(11)_z_loc[{m_},{j}]"
        )

# (12) RTAM_{k v} = sum_{i∈NF ∪ {o}} sum_{m∈MT} ET_m * y_{i m k v}
for k in K:
    for v in Vset:
        m.addConstr(
            RTAM[k,v] == quicksum(ET[m_]*y[i,m_,k,v] for i in (NF + [o]) for m_ in MT),
            name=f"(12)_RTAM[{k},{v}]"
        )

# (13) Time feasibility maintenance->flight: RTAM_{k v} <= DT_j if z_{m j k v} = 1
for m_ in MT:
    for j in NF:
        for k in K:
            for v in Vset:
                m.addConstr(
                    RTAM[k,v] - DT[j] <= Mbig * (1 - z[m_,j,k,v]),
                    name=f"(13)_time_z[{m_},{j},{k},{v}]"
                )

# (14) Max takeoffs between maint: sum_{i∈NF∪{o}} sum_{j∈NF} x_{i j k v} <= Cmax
for k in K:
    for v in Vset:
        m.addConstr(
            quicksum(x[i,j,k,v] for i in (NF + [o]) for j in NF) <= Cmax,
            name=f"(14)_takeoffs[{k},{v}]"
        )

# (15) Max flying time (first visit v=1): sum_{i∈NF∪{o}} sum_{j∈NF} DT_j * x_{i j k 1} <= Tmax
for k in K:
    m.addConstr(
        quicksum(DT[j] * x[i,j,k,1] for i in (NF + [o]) for j in NF) <= Tmax,
        name=f"(15)_time_first[{k}]"
    )

# (16) Max flying time for later visits v>1:
#      sum DT_j * x_{i j k v} + sum DT_j * z_{m j k v} <= Tmax
for k in K:
    for v in Vset:
        if v != 1:
            m.addConstr(
                quicksum(DT[j]*x[i,j,k,v] for i in (NF + [o]) for j in NF)
              + quicksum(DT[j]*z[m_,j,k,v] for m_ in MT for j in NF)
              <= Tmax,
              name=f"(16)_time_later[{k},{v}]"
            )

# (17) Each aircraft performs exactly Vmax maintenance operations
for k in K:
    m.addConstr(
        quicksum(y[i,m_,k,v] for i in (NF + [o]) for m_ in MT for v in Vset) == Vmax,
        name=f"(17)_maint_count[{k}]"
    )

# (18) V ≥ 1  — enforced by construction of Vset

# (19) Workforce capacity at each maintenance station (HARD CAP)
for m_ in MT:
    m.addConstr(
        quicksum(y[i,m_,k,v] for i in (NF + [o]) for k in K for v in Vset) <= MP[m_],
        name=f"(19)_workforce_cap[{m_}]"
    )

# (20)–(22) variable domains already declared as binaries/continuous

# (23) RTAM_{k v} > 0 — model strictness with a tiny epsilon if desired
eps = 1e-6
for k in K:
    for v in Vset:
        m.addConstr(RTAM[k,v] >= eps, name=f"(23)_RTAM_pos[{k},{v}]")

# =========================
# Optimize
# =========================
m.Params.OutputFlag = 1
m.optimize()

# =========================
# Extract solution (example)
# =========================
if m.status in (GRB.OPTIMAL, GRB.TIME_LIMIT):
    print("Objective value:", m.objVal)
    # Access values as x[i,j,k,v].X, etc.

