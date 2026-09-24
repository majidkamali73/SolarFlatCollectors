# -*- coding: utf-8 -*-
"""
Solar collector model v0.10.2D
Coupled engineering 2-D absorber plate + 1-D fluid model.
Python 3.7 compatible; standard library only.

This is an engineering model, not CFD. It resolves absorber-plate temperature
in x-y and couples each tube to its own fluid temperature. Header flow
maldistribution is retained from the v0.9 engineering hydraulic model.
"""
import math

SIGMA = 5.670374419e-8

# ---------------- hydraulic core ----------------
def friction_factor_laminar(Re):
    return 64.0 / Re if Re > 0.0 else 0.0

def friction_factor_haaland(Re, rel_roughness=0.0):
    if Re <= 0.0:
        return 0.0
    rr = max(rel_roughness, 0.0)
    return (-1.8 * math.log10((rr / 3.7) ** 1.11 + 6.9 / Re)) ** -2

def friction_factor(Re, rel_roughness=0.0):
    if Re <= 0.0:
        return 0.0, "zero_flow"
    if Re < 2300.0:
        return friction_factor_laminar(Re), "laminar"
    f2300 = friction_factor_laminar(2300.0)
    f4000 = friction_factor_haaland(4000.0, rel_roughness)
    if Re < 4000.0:
        x = (Re - 2300.0) / 1700.0
        return f2300 + x * (f4000 - f2300), "transition_2300_4000"
    return friction_factor_haaland(Re, rel_roughness), "turbulent_haaland"

def tube_dp(mdot, di, rho, mu, L, rel_roughness=0.0):
    if mdot <= 0.0:
        return 0.0, 0.0, 0.0, "zero_flow"
    A = math.pi * di ** 2 / 4.0
    V = mdot / (rho * A)
    Re = rho * V * di / mu
    f, regime = friction_factor(Re, rel_roughness)
    dp = f * (L / di) * 0.5 * rho * V ** 2
    return dp, V, Re, regime

def header_section_dp(mdot, D, rho, mu, length, rel_roughness=0.0):
    if mdot <= 0.0 or length <= 0.0:
        return 0.0
    A = math.pi * D ** 2 / 4.0
    V = mdot / (rho * A)
    Re = rho * V * D / mu
    f, _ = friction_factor(Re, rel_roughness)
    return f * (length / D) * 0.5 * rho * V ** 2

def _header_path_dp(branch_index, flows, connection_fraction, spacing,
                    Dh, rho, mu, rel_roughness=0.0):
    N = len(flows)
    if N <= 1:
        return 0.0
    x = [i * spacing for i in range(N)]
    xc = connection_fraction * (N - 1) * spacing
    xb = x[branch_index]
    if abs(xb - xc) < 1e-14:
        return 0.0
    if xb > xc:
        indices = [i for i in range(N) if x[i] >= xc and x[i] <= xb]
        indices.sort(key=lambda i: x[i])
        direction = 1
    else:
        indices = [i for i in range(N) if x[i] <= xc and x[i] >= xb]
        indices.sort(key=lambda i: x[i], reverse=True)
        direction = -1
    dp = 0.0
    for pos, j in enumerate(indices):
        xj = x[j]
        x0 = xc if pos == 0 else x[indices[pos - 1]]
        seg_len = abs(xj - x0)
        if direction > 0:
            remaining = sum(flows[k] for k in range(N) if x[k] >= xj - 1e-14)
        else:
            remaining = sum(flows[k] for k in range(N) if x[k] <= xj + 1e-14)
        if seg_len > 0.0 and remaining > 0.0:
            dp += header_section_dp(remaining, Dh, rho, mu, seg_len, rel_roughness)
        if j == branch_index:
            break
    return dp

def _path_values(flows, di, Dh, rho, mu, L_tube, spacing,
                 inlet_fraction, outlet_fraction, K_in, K_out,
                 rel_roughness):
    paths = []
    N = len(flows)
    for i in range(N):
        dp_hi = _header_path_dp(i, flows, inlet_fraction, spacing, Dh, rho, mu, rel_roughness)
        dp_ho = _header_path_dp(i, flows, outlet_fraction, spacing, Dh, rho, mu, rel_roughness)
        dp_t, V, Re, reg = tube_dp(flows[i], di, rho, mu, L_tube, rel_roughness)
        dp_b = (K_in + K_out) * 0.5 * rho * V ** 2
        paths.append((dp_hi + dp_t + dp_ho + dp_b, dp_hi, dp_ho, dp_t, dp_b, V, Re, reg))
    return paths

def solve_header_distribution(mdot_total, N, di, rho, mu, L_tube, spacing, Dh,
                              inlet_fraction=0.5, outlet_fraction=0.5,
                              K_branch_in=0.0, K_branch_out=0.0,
                              rel_roughness=0.0, max_iter=5000, tol=1e-7):
    target = mdot_total / float(N)
    q = [target] * N
    for _ in range(max_iter):
        paths = _path_values(q, di, Dh, rho, mu, L_tube, spacing,
                             inlet_fraction, outlet_fraction, K_branch_in, K_branch_out,
                             rel_roughness)
        dp = [p[0] for p in paths]
        ref = sum(dp) / N
        newq = [q[i] * math.sqrt(max(ref, 1e-12) / max(dp[i], 1e-12)) for i in range(N)]
        scale = mdot_total / sum(newq)
        newq = [x * scale for x in newq]
        err = max(abs(newq[i] - q[i]) for i in range(N)) / target
        q = newq
        if err < tol:
            break
    paths = _path_values(q, di, Dh, rho, mu, L_tube, spacing,
                         inlet_fraction, outlet_fraction, K_branch_in, K_branch_out,
                         rel_roughness)
    return {
        "tube_flows": q,
        "flow_nonuniformity": (max(q) - min(q)) / (sum(q) / N),
        "dp_path_Pa": sum(p[0] for p in paths) / N,
        "paths": paths
    }

# ---------------- thermal correlations ----------------
def nusselt_v06(Re, Pr):
    if Re < 2300.0:
        return 4.36, "laminar_constant_heat_flux"
    f4000 = friction_factor_haaland(4000.0, 0.0)
    nu4000 = ((f4000 / 8.0) * (4000.0 - 1000.0) * Pr /
              (1.0 + 12.7 * math.sqrt(f4000 / 8.0) * (Pr ** (2.0 / 3.0) - 1.0)))
    if Re < 4000.0:
        gamma = (Re - 2300.0) / 1700.0
        return 4.36 + gamma * (nu4000 - 4.36), "transition_2300_4000"
    f = friction_factor_haaland(Re, 0.0)
    nu = ((f / 8.0) * (Re - 1000.0) * Pr /
          (1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1.0)))
    return nu, "turbulent_gnielinski"

def top_loss_coefficient(Tpm_C, Ta_C, hw, M, C, f, eps_p, eps_c):
    Tpm = Tpm_C + 273.15
    Ta = Ta_C + 273.15
    dT = max(Tpm - Ta, 0.05)
    conv_res = M / ((C / Tpm) * ((dT / (M + f)) ** 0.33)) + 1.0 / hw
    Ut_conv = 1.0 / conv_res
    den_rad = (1.0 / (eps_p + 0.05 * M * (1.0 - eps_p)) +
               (2.0 * M + f - 1.0) / eps_c - M)
    Ut_rad = SIGMA * (Tpm ** 2 + Ta ** 2) * (Tpm + Ta) / den_rad
    return Ut_conv + Ut_rad

def thermal_tube(mdot, di, rho, mu, k, cp):
    A = math.pi * di ** 2 / 4.0
    V = mdot / (rho * A)
    Re = rho * V * di / mu
    Pr = cp * mu / k
    Nu, regime = nusselt_v06(Re, Pr)
    hf = Nu * k / di
    return V, Re, Pr, Nu, hf, regime

# ---------------- coupled model ----------------
def coupled_case(mdot_total=0.1, N=10, w=0.10, L1=3.0,
                 di=0.370*0.0254, do=0.500*0.0254,
                 plate_k=237.0, plate_delta=0.0005,
                 adhesive_k=0.2, adhesive_delta=0.0005,
                 insulation_k=0.04, deltab=0.04,
                 solar_I=760.0, Ta=15.0, Tfi=20.0, wind=2.0,
                 beta=30.0, M=1.0, tau=0.90, n_cover=1.52,
                 eps_p=0.95, alpha=0.95, Dh_ratio=3.0,
                 connection_fraction=0.5, K_branch_in=0.0, K_branch_out=0.0,
                 pump_efficiency=0.65, nx=36, ny_per_gap=5,
                 max_outer=35, max_plate=800, relax=0.55):
    rho, mu, cp, k = 998.0, 0.0004275, 4180.0, 0.60
    L2 = N * w
    Ap = L1 * L2
    dx = L1 / float(nx)
    ny = max(N * ny_per_gap, N + 2)
    dy = L2 / float(ny)
    S = solar_I * (tau ** M) * alpha
    hw = 5.7 + 3.8 * wind
    fwind = (1.0 - 0.04 * hw + 0.0005 * hw ** 2) * (1.0 + 0.091 * M)
    C = 365.9 * (1.0 - 0.00883 * beta + 0.0001298 * beta ** 2)
    R = ((1.0 - n_cover) / (2.0 + n_cover)) ** 2
    eps_c = 1.0 - R - tau
    Ub = insulation_k / deltab
    L3 = deltab + M * 0.03 + 0.01
    Us = ((L1 + L2) * L3 * insulation_k) / (L1 * L2 * deltab)

    hyd = solve_header_distribution(mdot_total, N, di, rho, mu, L1, w,
                                    Dh_ratio * di, connection_fraction,
                                    connection_fraction, K_branch_in, K_branch_out)
    flows = hyd["tube_flows"]
    tube_data = [thermal_tube(m, di, rho, mu, k, cp) for m in flows]

    # Tube center positions; grid rows are assigned to the nearest center.
    ycenters = [(i + 0.5) * w for i in range(N)]
    tube_rows = [min(range(ny), key=lambda j: abs((j + 0.5) * dy - yc)) for yc in ycenters]
    tube_row_map = {}
    for i, row in enumerate(tube_rows):
        tube_row_map.setdefault(row, []).append(i)

    # Initial plate temperature and fluid profile.
    T = [[Tfi for _ in range(ny)] for _ in range(nx)]
    tube_tf = [[Tfi for _ in range(nx + 1)] for _ in range(N)]
    UL = 4.0

    for outer in range(max_outer):
        old_T = [r[:] for r in T]
        Tmean = sum(sum(r) for r in T) / (nx * ny)
        UL = top_loss_coefficient(Tmean, Ta, hw, M, C, fwind, eps_p, eps_c) + Ub + Us

        # March each tube with the current plate temperature.
        for i in range(N):
            Gprime = 1.0 / (1.0 / (tube_data[i][4] * math.pi * di) +
                             adhesive_delta / (adhesive_k * math.pi * do))
            tube_tf[i][0] = Tfi
            mcp = flows[i] * cp
            for ix in range(nx):
                Tp = T[ix][tube_rows[i]]
                a = Gprime * dx / mcp
                e = math.exp(-a)
                tube_tf[i][ix + 1] = Tp - (Tp - tube_tf[i][ix]) * e

        # Gauss-Seidel plate solve with the fluid temperatures held fixed.
        Gx = plate_k * plate_delta * dy / dx
        Gy = plate_k * plate_delta * dx / dy
        cell_area = dx * dy
        for _ in range(max_plate):
            max_change = 0.0
            for ix in range(nx):
                for iy in range(ny):
                    diag = UL * cell_area
                    rhs = S * cell_area + UL * cell_area * Ta
                    if ix > 0:
                        diag += Gx; rhs += Gx * T[ix - 1][iy]
                    if ix < nx - 1:
                        diag += Gx; rhs += Gx * T[ix + 1][iy]
                    if iy > 0:
                        diag += Gy; rhs += Gy * T[ix][iy - 1]
                    if iy < ny - 1:
                        diag += Gy; rhs += Gy * T[ix][iy + 1]
                    if iy in tube_row_map:
                        for ti in tube_row_map[iy]:
                            Gprime = 1.0 / (1.0 / (tube_data[ti][4] * math.pi * di) +
                                             adhesive_delta / (adhesive_k * math.pi * do))
                            gcell = Gprime * dx
                            diag += gcell
                            rhs += gcell * tube_tf[ti][ix]
                    newv = rhs / diag
                    dv = newv - T[ix][iy]
                    T[ix][iy] += relax * dv
                    if abs(dv) > max_change:
                        max_change = abs(dv)
            if max_change < 1e-6:
                break

        diff = max(abs(T[ix][iy] - old_T[ix][iy]) for ix in range(nx) for iy in range(ny))
        if diff < 2e-5:
            break

    # Final fluid march and heat output.
    q_tubes = []
    Tfo = []
    for i in range(N):
        Gprime = 1.0 / (1.0 / (tube_data[i][4] * math.pi * di) +
                         adhesive_delta / (adhesive_k * math.pi * do))
        tube_tf[i][0] = Tfi
        mcp = flows[i] * cp
        for ix in range(nx):
            Tp = T[ix][tube_rows[i]]
            a = Gprime * dx / mcp
            tube_tf[i][ix + 1] = Tp - (Tp - tube_tf[i][ix]) * math.exp(-a)
        Tfo.append(tube_tf[i][-1])
        q_tubes.append(flows[i] * cp * (Tfo[-1] - Tfi))
    q = sum(q_tubes)
    eta = q / (Ap * solar_I)
    Tavg = sum(sum(r) for r in T) / (nx * ny)
    Tmax = max(max(r) for r in T)
    Tmin = min(min(r) for r in T)
    Qvol = mdot_total / rho
    pump_power = hyd["dp_path_Pa"] * Qvol / pump_efficiency
    absorbed = S * Ap
    losses = sum(sum(UL * (T[ix][iy] - Ta) * cell_area for iy in range(ny)) for ix in range(nx))
    balance_error = (absorbed - losses - q) / max(abs(absorbed), 1.0)
    return {
        "N": N, "w_m": w, "L1_m": L1, "L2_m": L2, "Ap_m2": Ap,
        "q_W": q, "efficiency": eta, "UL_W_m2K": UL,
        "Tavg_C": Tavg, "Tmax_C": Tmax, "Tmin_C": Tmin,
        "Tfo_mean_C": sum(Tfo) / N, "Tfo_min_C": min(Tfo), "Tfo_max_C": max(Tfo),
        "Re_mean": sum(x[1] for x in tube_data) / N,
        "hf_mean": sum(x[4] for x in tube_data) / N,
        "flow_nonuniformity": hyd["flow_nonuniformity"],
        "dp_system_Pa": hyd["dp_path_Pa"], "pump_power_W": pump_power,
        "balance_error_fraction": balance_error,
        "tube_flows": flows, "tube_Tfo": Tfo,
        "plate_temperature_C": T,
        "iterations_outer": outer + 1
    }

def solve_L1_for_q(mdot_total, N, w, q_target, L1_lo=1.0, L1_hi=5.0, **kwargs):
    # Expand upper bound if necessary.
    rlo = coupled_case(mdot_total, N, w, L1_lo, **kwargs)
    rhi = coupled_case(mdot_total, N, w, L1_hi, **kwargs)
    for _ in range(6):
        if rhi["q_W"] >= q_target:
            break
        L1_hi *= 1.4
        rhi = coupled_case(mdot_total, N, w, L1_hi, **kwargs)
    if rlo["q_W"] >= q_target:
        return rlo
    if rhi["q_W"] < q_target:
        raise ValueError("Target duty not reached by L1 upper bound")
    for _ in range(18):
        mid = 0.5 * (L1_lo + L1_hi)
        rm = coupled_case(mdot_total, N, w, mid, **kwargs)
        if rm["q_W"] < q_target:
            L1_lo = mid
        else:
            L1_hi = mid
    return coupled_case(mdot_total, N, w, 0.5 * (L1_lo + L1_hi), **kwargs)


def coupled_case_fast(mdot_total=0.1, N=10, w=0.10, L1=3.0,
                 di=0.370*0.0254, do=0.500*0.0254,
                 plate_k=237.0, plate_delta=0.0005,
                 adhesive_k=0.2, adhesive_delta=0.0005,
                 insulation_k=0.04, deltab=0.04,
                 solar_I=760.0, Ta=15.0, Tfi=20.0, wind=2.0,
                 beta=30.0, M=1.0, tau=0.90, n_cover=1.52,
                 eps_p=0.95, alpha=0.95, Dh_ratio=3.0,
                 connection_fraction=0.5, K_branch_in=0.0, K_branch_out=0.0,
                 pump_efficiency=0.65, nx=40, ny_per_gap=5,
                 max_outer=30, relax=0.6):
    # Pure-Python implementation: intentionally no NumPy/SciPy dependency.
    # This is compatible with stock Python 3.7.2 + IDLE.
    rho, mu, cp, k = 998.0, 0.0004275, 4180.0, 0.60
    L2=N*w; Ap=L1*L2; dx=L1/float(nx); ny=max(N*ny_per_gap,N+2); dy=L2/float(ny)
    S=solar_I*(tau**M)*alpha; hw=5.7+3.8*wind
    fwind=(1-.04*hw+.0005*hw**2)*(1+.091*M); C=365.9*(1-.00883*beta+.0001298*beta**2)
    R=((1-n_cover)/(2+n_cover))**2; eps_c=1-R-tau
    Ub=insulation_k/deltab; L3=deltab+M*.03+.01
    Us=((L1+L2)*L3*insulation_k)/(L1*L2*deltab)
    hyd=solve_header_distribution(mdot_total,N,di,rho,mu,L1,w,Dh_ratio*di,connection_fraction,connection_fraction,K_branch_in,K_branch_out)
    flows=hyd['tube_flows']
    tube_data=[thermal_tube(m,di,rho,mu,k,cp) for m in flows]
    tube_rows=[min(range(ny),key=lambda j:abs((j+.5)*dy-(i+.5)*w)) for i in range(N)]
    T=[[Tfi for _ in range(ny)] for _ in range(nx)]
    tube_tf=[[Tfi for _ in range(nx+1)] for _ in range(N)]
    cell_area=dx*dy
    UL=4.0
    for outer in range(max_outer):
        max_change=0.0
        Tprev=[row[:] for row in T]
        Tmean=sum(sum(row) for row in T)/(nx*ny)
        UL=top_loss_coefficient(Tmean,Ta,hw,M,C,fwind,eps_p,eps_c)+Ub+Us
        # Fluid march using the previous/current plate temperatures.
        for i in range(N):
            hf=tube_data[i][4]
            Gp=1/(1/(hf*math.pi*di)+adhesive_delta/(adhesive_k*math.pi*do))
            mcp=flows[i]*cp
            a=Gp*dx/mcp; e=math.exp(-a); tube_tf[i][0]=Tfi
            for ix in range(nx):
                tube_tf[i][ix+1]=T[ix][tube_rows[i]]-(T[ix][tube_rows[i]]-tube_tf[i][ix])*e
        # Gauss-Seidel solution of the 2-D plate conduction equation.
        Gx=plate_k*plate_delta*dy/dx
        Gy=plate_k*plate_delta*dx/dy
        for sweep in range(350):
            sweep_change=0.0
            for ix in range(nx):
                for iy in range(ny):
                    diag=UL*cell_area
                    rhs=S*cell_area+UL*cell_area*Ta
                    if ix>0:
                        diag+=Gx; rhs+=Gx*T[ix-1][iy]
                    if ix<nx-1:
                        diag+=Gx; rhs+=Gx*T[ix+1][iy]
                    if iy>0:
                        diag+=Gy; rhs+=Gy*T[ix][iy-1]
                    if iy<ny-1:
                        diag+=Gy; rhs+=Gy*T[ix][iy+1]
                    for ti,row in enumerate(tube_rows):
                        if row==iy:
                            hf=tube_data[ti][4]
                            Gp=1/(1/(hf*math.pi*di)+adhesive_delta/(adhesive_k*math.pi*do))
                            a=Gp*dx/(flows[ti]*cp)
                            gc=flows[ti]*cp*(1-math.exp(-a))
                            diag+=gc; rhs+=gc*tube_tf[ti][ix]
                    newT=rhs/diag
                    d=abs(newT-T[ix][iy])
                    if d>sweep_change: sweep_change=d
                    T[ix][iy]=newT
            if sweep_change<2e-7: break
        for ix in range(nx):
            for iy in range(ny):
                d=abs(T[ix][iy]-Tprev[ix][iy])
                if d>max_change: max_change=d
        if max_change<2e-5: break
    # final fluid march and outputs
    q_t=[]; Tfo=[]
    for i in range(N):
        hf=tube_data[i][4]; Gp=1/(1/(hf*math.pi*di)+adhesive_delta/(adhesive_k*math.pi*do)); mcp=flows[i]*cp; a=Gp*dx/mcp; e=math.exp(-a); tf=Tfi
        for ix in range(nx): tf=T[ix][tube_rows[i]]-(T[ix][tube_rows[i]]-tf)*e
        Tfo.append(tf); q_t.append(flows[i]*cp*(tf-Tfi))
    q=sum(q_t); eta=q/(Ap*solar_I); Tavg=sum(sum(row) for row in T)/(nx*ny); Tmax=max(max(row) for row in T); Tmin=min(min(row) for row in T)
    pump=hyd['dp_path_Pa']*(mdot_total/rho)/pump_efficiency
    absorbed=S*Ap; losses=sum(UL*(T[ix][iy]-Ta)*cell_area for ix in range(nx) for iy in range(ny)); bal=(absorbed-losses-q)/max(abs(absorbed),1.)
    return {'N':N,'w_m':w,'L1_m':L1,'L2_m':L2,'Ap_m2':Ap,'q_W':q,'efficiency':eta,'UL_W_m2K':UL,'Tavg_C':Tavg,'Tmax_C':Tmax,'Tmin_C':Tmin,'Tfo_mean_C':sum(Tfo)/N,'Tfo_min_C':min(Tfo),'Tfo_max_C':max(Tfo),'Re_mean':sum(x[1] for x in tube_data)/N,'hf_mean':sum(x[4] for x in tube_data)/N,'flow_nonuniformity':hyd['flow_nonuniformity'],'dp_system_Pa':hyd['dp_path_Pa'],'pump_power_W':pump,'balance_error_fraction':bal,'tube_flows':flows,'tube_Tfo':Tfo,'plate_temperature_C':T,'iterations_outer':outer+1}

def solve_L1_for_q_fast(mdot_total,N,w,q_target,L1_lo=1.0,L1_hi=5.0,**kwargs):
    rlo=coupled_case_fast(mdot_total,N,w,L1_lo,**kwargs); rhi=coupled_case_fast(mdot_total,N,w,L1_hi,**kwargs)
    while rhi['q_W']<q_target:
        L1_hi*=1.4; rhi=coupled_case_fast(mdot_total,N,w,L1_hi,**kwargs)
        if L1_hi>15: raise ValueError('Target duty not reached')
    if rlo['q_W']>=q_target: return rlo
    for _ in range(14):
        mid=(L1_lo+L1_hi)/2; rm=coupled_case_fast(mdot_total,N,w,mid,**kwargs)
        if rm['q_W']<q_target: L1_lo=mid
        else: L1_hi=mid
    return coupled_case_fast(mdot_total,N,w,(L1_lo+L1_hi)/2,**kwargs)
