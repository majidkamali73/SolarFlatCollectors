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

def _header_dps(flows, connection_fraction, spacing, Dh, rho, mu,
                rel_roughness=0.0):
    # Header pressure loss from the connection point to every branch.
    # Each header segment carries the flow of all branches beyond it.
    N = len(flows)
    out = [0.0] * N
    if N <= 1:
        return out
    xc = connection_fraction * (N - 1) * spacing
    beyond = [0.0] * N
    total = 0.0
    for j in range(N - 1, -1, -1):
        total += flows[j]
        beyond[j] = total
    dp = 0.0
    prev = xc
    for j in range(N):
        xj = j * spacing
        if xj < xc:
            continue
        if xj - prev > 0.0 and beyond[j] > 0.0:
            dp += header_section_dp(beyond[j], Dh, rho, mu, xj - prev, rel_roughness)
        prev = xj
        if xj - xc >= 1e-14:
            out[j] = dp
    total = 0.0
    for j in range(N):
        total += flows[j]
        beyond[j] = total
    dp = 0.0
    prev = xc
    for j in range(N - 1, -1, -1):
        xj = j * spacing
        if xj > xc:
            continue
        if prev - xj > 0.0 and beyond[j] > 0.0:
            dp += header_section_dp(beyond[j], Dh, rho, mu, prev - xj, rel_roughness)
        prev = xj
        if xc - xj >= 1e-14:
            out[j] = dp
    return out

def _path_values(flows, di, Dh, rho, mu, L_tube, spacing,
                 inlet_fraction, outlet_fraction, K_in, K_out,
                 rel_roughness):
    paths = []
    N = len(flows)
    dp_in = _header_dps(flows, inlet_fraction, spacing, Dh, rho, mu, rel_roughness)
    dp_out = _header_dps(flows, outlet_fraction, spacing, Dh, rho, mu, rel_roughness)
    for i in range(N):
        dp_hi = dp_in[i]
        dp_ho = dp_out[i]
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

# ---------------- water properties ----------------
# Saturated liquid water, 0-100 C, in 5 K steps (standard property-table
# values): density kg/m3, specific heat J/kg.K, conductivity W/m.K,
# dynamic viscosity Pa.s. Linear interpolation; clamped outside the range.
_WATER_TABLE = (
    (0.0, 999.8, 4217.0, 0.561, 1.792e-3),
    (5.0, 999.9, 4205.0, 0.571, 1.519e-3),
    (10.0, 999.7, 4194.0, 0.580, 1.307e-3),
    (15.0, 999.1, 4185.0, 0.589, 1.138e-3),
    (20.0, 998.0, 4182.0, 0.598, 1.002e-3),
    (25.0, 997.0, 4180.0, 0.607, 0.891e-3),
    (30.0, 996.0, 4178.0, 0.615, 0.798e-3),
    (35.0, 994.0, 4178.0, 0.623, 0.720e-3),
    (40.0, 992.1, 4179.0, 0.631, 0.653e-3),
    (45.0, 990.1, 4180.0, 0.637, 0.596e-3),
    (50.0, 988.1, 4181.0, 0.644, 0.547e-3),
    (55.0, 985.2, 4183.0, 0.649, 0.504e-3),
    (60.0, 983.3, 4185.0, 0.654, 0.467e-3),
    (65.0, 980.4, 4187.0, 0.659, 0.433e-3),
    (70.0, 977.5, 4190.0, 0.663, 0.404e-3),
    (75.0, 974.7, 4193.0, 0.667, 0.378e-3),
    (80.0, 971.8, 4197.0, 0.670, 0.355e-3),
    (85.0, 968.1, 4201.0, 0.673, 0.333e-3),
    (90.0, 965.3, 4206.0, 0.675, 0.315e-3),
    (95.0, 961.5, 4212.0, 0.677, 0.297e-3),
    (100.0, 957.9, 4217.0, 0.679, 0.282e-3),
)

def water_properties(T_C):
    """Return (rho, mu, cp, k) of liquid water at T_C."""
    T = min(max(T_C, _WATER_TABLE[0][0]), _WATER_TABLE[-1][0])
    for lo, hi in zip(_WATER_TABLE[:-1], _WATER_TABLE[1:]):
        if T <= hi[0]:
            x = (T - lo[0]) / (hi[0] - lo[0])
            rho, cp, k, mu = [lo[j] + x * (hi[j] - lo[j]) for j in range(1, 5)]
            return rho, mu, cp, k

def _fluid_state(Tf_C, mdot_total, N, di, L1, w, Dh, connection_fraction,
                 K_branch_in, K_branch_out):
    # Properties, header flow split and tube-side coefficients at Tf_C.
    rho, mu, cp, k = water_properties(Tf_C)
    hyd = solve_header_distribution(mdot_total, N, di, rho, mu, L1, w, Dh,
                                    connection_fraction, connection_fraction,
                                    K_branch_in, K_branch_out)
    flows = hyd["tube_flows"]
    tube_data = [thermal_tube(m, di, rho, mu, k, cp) for m in flows]
    return rho, mu, cp, k, hyd, flows, tube_data

# ---------------- coupled model ----------------
CONNECTION_TYPES = ('below_plate', 'above_plate', 'in_line')


def coupled_case(*args, **kwargs):
    # Kept for older scripts; the plate solver options of the former
    # point-iteration version are accepted and ignored.
    kwargs.pop('max_plate', None)
    return coupled_case_fast(*args, **kwargs)


def solve_L1_for_q(mdot_total, N, w, q_target, L1_lo=1.0, L1_hi=5.0, **kwargs):
    kwargs.pop('max_plate', None)
    return solve_L1_for_q_fast(mdot_total, N, w, q_target, L1_lo, L1_hi, **kwargs)


def coupled_case_fast(mdot_total=0.1, N=10, w=0.10, L1=3.0,
                 di=0.370*0.0254, do=0.500*0.0254,
                 plate_k=237.0, plate_delta=0.0005,
                 adhesive_k=0.2, adhesive_delta=0.0005,
                 insulation_k=0.04, deltab=0.04,
                 side_insulation_k=None, side_insulation_thickness_m=None,
                 solar_I=760.0, Ta=15.0, Tfi=20.0, wind=2.0,
                 beta=30.0, M=1.0, tau=0.90, n_cover=1.52,
                 eps_p=0.95, alpha=0.95, Dh_ratio=3.0,
                 connection_fraction=0.5, K_branch_in=0.0, K_branch_out=0.0,
                 pump_efficiency=0.65, nx=40, ny_per_gap=5,
                 max_outer=30, relax=0.6, cover_emissivity=None,
                 connection_type='below_plate'):
    # Pure-Python implementation: intentionally no NumPy/SciPy dependency.
    # This is compatible with stock Python 3.7.2 + IDLE.
    #
    # Tube-to-plate connection, following the three legacy arrangements:
    #   below_plate: plate -> bond -> tube wall -> water; the bond is as wide
    #                as the tube, so its resistance per metre of tube is
    #                adhesive_delta / (adhesive_k * do).
    #   above_plate: the strip of width do under the tube gains and loses
    #                heat at the tube-wall temperature; only the heat coming
    #                from the fins crosses the bond.
    #   in_line:     the tube is part of the plate, no bond resistance.
    # In all three the plate strip of width do at the tube carries no fin
    # resistance, as in the classical F-prime relations.
    if connection_type not in CONNECTION_TYPES:
        raise ValueError('connection_type must be below_plate, above_plate or in_line')
    L2=N*w; Ap=L1*L2; dx=L1/float(nx); ny=max(N*ny_per_gap,N+2); dy=L2/float(ny)
    S=solar_I*(tau**M)*alpha; hw=5.7+3.8*wind
    fwind=(1-.04*hw+.0005*hw**2)*(1+.091*M); C=365.9*(1-.00883*beta+.0001298*beta**2)
    eps_c=0.88 if cover_emissivity is None else cover_emissivity
    Ub=insulation_k/deltab
    side_k=insulation_k if side_insulation_k is None else side_insulation_k
    side_delta=deltab if side_insulation_thickness_m is None else side_insulation_thickness_m
    if side_k<=0 or side_delta<=0:
        raise ValueError('Side insulation conductivity and thickness must be positive.')
    L3=deltab+M*.03+.01
    Us=((L1+L2)*L3*side_k)/(L1*L2*side_delta)
    Tf_prop=Tfi
    rho,mu,cp,k,hyd,flows,tube_data=_fluid_state(Tf_prop,mdot_total,N,di,L1,w,Dh_ratio*di,connection_fraction,K_branch_in,K_branch_out)
    tube_rows=[min(range(ny),key=lambda j:abs((j+.5)*dy-(i+.5)*w)) for i in range(N)]
    is_tube_row=[False]*ny
    for row in tube_rows: is_tube_row[row]=True
    T=[[Tfi for _ in range(ny)] for _ in range(nx)]
    tube_tf=[[Tfi for _ in range(nx+1)] for _ in range(N)]
    cell_area=dx*dy
    strip=min(do,dy)
    above=connection_type=='above_plate'
    bond_R=0.0 if connection_type=='in_line' else adhesive_delta/(adhesive_k*do)
    # Conductance of each y-face; next to a tube the conduction path is
    # shorter by half the strip width.
    Gx=plate_k*plate_delta*dy/dx
    gy=[0.0]*ny
    for iy in range(ny-1):
        cut=0.5*strip*(is_tube_row[iy]+is_tube_row[iy+1])
        gy[iy]=plate_k*plate_delta*dx/max(dy-cut,0.25*dy)
    UL=4.0

    def tube_terms(ti):
        # Returns (a, lam, G1, G2, G3) of one tube: the water approaches
        # lam*Tplate + (1-lam)*Tstar with exponent a per x-step.
        G2=tube_data[ti][4]*math.pi*di
        mcp=flows[ti]*cp
        if above:
            G1=1.0/bond_R; G3=UL*strip
            return G2*(G1+G3)/(G1+G2+G3)*dx/mcp, G1/(G1+G3), G1, G2, G3
        return dx/((1.0/G2+bond_R)*mcp), 1.0, 0.0, G2, 0.0

    def march():
        Tstar=Ta+S/UL
        for i in range(N):
            a,lam,G1,G2,G3=tube_terms(i); e=math.exp(-a); row=tube_rows[i]
            tf=Tfi; tube_tf[i][0]=Tfi
            for ix in range(nx):
                teq=lam*T[ix][row]+(1.0-lam)*Tstar
                tf=teq-(teq-tf)*e
                tube_tf[i][ix+1]=tf

    for outer in range(max_outer):
        max_change=0.0
        Tprev=[row[:] for row in T]
        Tmean=sum(sum(row) for row in T)/(nx*ny)
        UL=top_loss_coefficient(Tmean,Ta,hw,M,C,fwind,eps_p,eps_c)+Ub+Us
        Tstar=Ta+S/UL
        # Fluid march using the previous/current plate temperatures.
        march()
        # Water properties follow the mean fluid temperature.
        Tf_mean=0.5*(Tfi+sum(flows[i]*tube_tf[i][-1] for i in range(N))/mdot_total)
        props_updated=abs(Tf_mean-Tf_prop)>0.25
        if props_updated:
            Tf_prop=Tf_mean
            rho,mu,cp,k,hyd,flows,tube_data=_fluid_state(Tf_prop,mdot_total,N,di,L1,w,Dh_ratio*di,connection_fraction,K_branch_in,K_branch_out)
        # Plate conduction: each x-station is a tridiagonal system in y and is
        # solved exactly; the stations are swept until the x-coupling settles.
        area=[cell_area]*ny
        gc_row=[0.0]*ny; gt=[[0.0]*ny for _ in range(nx)]
        for ti,row in enumerate(tube_rows):
            a,lam,G1,G2,G3=tube_terms(ti); e=math.exp(-a); mcp=flows[ti]*cp
            if above:
                c=(1.0-e)/a; sg=G1+G2+G3
                area[row]-=strip*dx
                gc_row[row]+=dx*G1/sg*(G2+G3-G2*(1.0-c)*lam)
                kstar=dx*G1/sg*(G3+G2*(1.0-c)*(1.0-lam))*Tstar; kf=dx*G1/sg*G2*c
                for ix in range(nx): gt[ix][row]+=kstar+kf*tube_tf[ti][ix]
            else:
                gc=mcp*(1.0-e)
                gc_row[row]+=gc
                for ix in range(nx): gt[ix][row]+=gc*tube_tf[ti][ix]
        # Elimination factors for a station with 0, 1 or 2 x-neighbours.
        factors=[]
        for n_side in (0,1,2):
            inv=[0.0]*ny; up=[0.0]*ny; prev_up=0.0
            for iy in range(ny):
                diag=UL*area[iy]+n_side*Gx+gc_row[iy]+gy[iy]
                if iy>0: diag+=gy[iy-1]; diag-=gy[iy-1]*prev_up
                inv[iy]=1.0/diag
                prev_up=gy[iy]*inv[iy]; up[iy]=prev_up
            factors.append((inv,up))
        fwd=[0.0]*ny
        src=[(S+UL*Ta)*area[iy] for iy in range(ny)]
        for sweep in range(350):
            sweep_change=0.0
            for ix in range(nx):
                left=T[ix-1] if ix>0 else None
                right=T[ix+1] if ix<nx-1 else None
                inv,up=factors[(left is not None)+(right is not None)]
                g=gt[ix]; acc=0.0
                for iy in range(ny):
                    rhs=src[iy]+g[iy]
                    if left is not None: rhs+=Gx*left[iy]
                    if right is not None: rhs+=Gx*right[iy]
                    if iy>0: rhs+=gy[iy-1]*acc
                    acc=rhs*inv[iy]
                    fwd[iy]=acc
                row=T[ix]; newT=0.0
                for iy in range(ny-1,-1,-1):
                    newT=fwd[iy]+up[iy]*newT if iy<ny-1 else fwd[iy]
                    d=abs(newT-row[iy])
                    if d>sweep_change: sweep_change=d
                    row[iy]=newT
            if sweep_change<2e-7: break
        for ix in range(nx):
            for iy in range(ny):
                d=abs(T[ix][iy]-Tprev[ix][iy])
                if d>max_change: max_change=d
        if max_change<2e-5 and not props_updated: break
    # final fluid march and outputs
    march()
    Tfo=[tube_tf[i][-1] for i in range(N)]
    q_t=[flows[i]*cp*(Tfo[i]-Tfi) for i in range(N)]
    q=sum(q_t); eta=q/(Ap*solar_I); Tavg=sum(sum(row) for row in T)/(nx*ny); Tmax=max(max(row) for row in T); Tmin=min(min(row) for row in T)
    pump=hyd['dp_path_Pa']*(mdot_total/rho)/pump_efficiency
    area=[cell_area]*ny
    losses=0.0
    if above:
        Tstar=Ta+S/UL
        for ti,row in enumerate(tube_rows):
            a,lam,G1,G2,G3=tube_terms(ti); c=(1.0-math.exp(-a))/a
            area[row]-=strip*dx
            for ix in range(nx):
                teq=lam*T[ix][row]+(1.0-lam)*Tstar
                tf_avg=teq-(teq-tube_tf[ti][ix])*c
                tw=(G1*T[ix][row]+G2*tf_avg+G3*Tstar)/(G1+G2+G3)
                losses+=UL*strip*dx*(tw-Ta)
    absorbed=S*Ap; losses+=sum(UL*(T[ix][iy]-Ta)*area[iy] for ix in range(nx) for iy in range(ny)); bal=(absorbed-losses-q)/max(abs(absorbed),1.)
    return {'N':N,'w_m':w,'L1_m':L1,'L2_m':L2,'Ap_m2':Ap,'q_W':q,'efficiency':eta,'UL_W_m2K':UL,'Tavg_C':Tavg,'Tmax_C':Tmax,'Tmin_C':Tmin,'Tfo_mean_C':sum(Tfo)/N,'Tfo_min_C':min(Tfo),'Tfo_max_C':max(Tfo),'Re_mean':sum(x[1] for x in tube_data)/N,'hf_mean':sum(x[4] for x in tube_data)/N,'flow_nonuniformity':hyd['flow_nonuniformity'],'dp_system_Pa':hyd['dp_path_Pa'],'pump_power_W':pump,'balance_error_fraction':bal,'tube_flows':flows,'tube_Tfo':Tfo,'plate_temperature_C':T,'iterations_outer':outer+1,'water_temperature_C':Tf_prop,'water_density_kg_m3':rho,'water_viscosity_Pa_s':mu,'water_cp_J_kgK':cp,'water_k_W_mK':k,'connection_type':connection_type}

def solve_L1_for_q_fast(mdot_total,N,w,q_target,L1_lo=1.0,L1_hi=5.0,L1_limit=None,q_tol=0.1,**kwargs):
    # Shortest length that delivers q_target, never shorter than L1_lo.
    # L1_limit, when given, is the longest acceptable collector: the search
    # stops at once if even that length cannot deliver the duty.
    rlo=coupled_case_fast(mdot_total,N,w,L1_lo,**kwargs)
    if rlo['q_W']>=q_target: return rlo
    if L1_limit is not None:
        if L1_limit<=L1_lo: raise ValueError('Target duty not reached')
        L1_hi=L1_limit; rhi=coupled_case_fast(mdot_total,N,w,L1_hi,**kwargs)
        if rhi['q_W']<q_target: raise ValueError('Target duty not reached')
    else:
        rhi=coupled_case_fast(mdot_total,N,w,L1_hi,**kwargs)
        while rhi['q_W']<q_target:
            L1_hi*=1.4; rhi=coupled_case_fast(mdot_total,N,w,L1_hi,**kwargs)
            if L1_hi>15: raise ValueError('Target duty not reached')
    # False position with the Illinois correction; rhi always meets the duty.
    flo=rlo['q_W']-q_target; fhi=rhi['q_W']-q_target; side=0
    for _ in range(40):
        if rhi['q_W']-q_target<=q_tol or L1_hi-L1_lo<=1e-6: break
        mid=L1_hi-fhi*(L1_hi-L1_lo)/(fhi-flo)
        rm=coupled_case_fast(mdot_total,N,w,mid,**kwargs); fm=rm['q_W']-q_target
        if fm<0:
            L1_lo=mid; flo=fm
            if side<0: fhi*=0.5
            side=-1
        else:
            L1_hi=mid; fhi=fm; rhi=rm
            if side>0: flo*=0.5
            side=1
    return rhi
