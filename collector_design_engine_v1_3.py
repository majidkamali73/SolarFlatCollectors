# Collector Design Engine v1.5
# Python 3.7.2 compatible; standard library only.
# v1.2 adds multi-material catalogs, rear insulation, 1/2 covers,
# cover spacing, adhesive catalog, and the three legacy tube/plate connection choices.

import csv, os, sys, math
HERE=os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0,HERE)
from solar_collector_model_v0_10_2d import solve_L1_for_q_fast


def read_catalog(name):
    path=os.path.join(HERE,'data',name+'.csv')
    with open(path,'r',newline='') as f:
        return list(csv.DictReader(f))


def catalog(name):
    return read_catalog(name)


def get_by_id(rows,item_id):
    for r in rows:
        if r.get('id')==item_id:
            return r
    raise ValueError('Item not found: '+str(item_id))


def _f(r,key):
    return float(r[key])


def validate_tube(tube):
    od=_f(tube,'od_m'); di=_f(tube,'id_m')
    if od<=di or di<=0:
        raise ValueError('Tube geometry invalid: OD must be greater than ID.')
    return 0.5*(od-di)


def cover_set(cover1_id, cover2_id=None):
    rows=catalog('covers')
    c1=get_by_id(rows,cover1_id)
    covers=[c1]
    if cover2_id:
        covers.append(get_by_id(rows,cover2_id))
    return covers


def effective_cover_properties(covers):
    # One cover is the legacy-equivalent path. For two different covers,
    # this is an explicit engineering approximation: optical transmission
    # is multiplied and refractive index is averaged. The approximation is
    # reported so it is not mistaken for a rigorous multi-layer radiation model.
    tau=1.0
    for c in covers:
        tau *= _f(c,'transmittance')
    n=sum(_f(c,'refractive_index') for c in covers)/float(len(covers))
    eps=sum(_f(c,'emissivity') for c in covers)/float(len(covers))
    return tau,n,eps


def design_collector(heat_required,inlet_temperature,ambient_temperature,
                     solar_irradiance,wind_speed,beta=30.0,fluid='water',
                     absorber_id='al_0p5',tube_id='cu_0p5_bwg16',
                     cover1_id='lowiron_4mm',cover2_id=None,cover_gap_m=0.012,
                     back_insulation_id='glass_wool',back_insulation_thickness_m=None,
                     side_insulation_id='glass_wool',side_insulation_thickness_m=None,
                     insulation_id=None,insulation_thickness_m=None,
                     adhesive_id='default_ka',adhesive_k_W_mK=None,adhesive_thickness_m=None,
                     connection_type='below_plate',mdot_total=0.10,
                     N_values=range(8,21),spacing_values=(0.06,0.08,0.10,0.12,0.14,0.16),
                     max_L1=3.5,max_L2=1.5,max_nonuniformity_percent=5.0,
                     objective='min_cost',tube_price_per_m=None,
                     absorber_price_per_m2=None,cover_price_per_m2=None,
                     header_price_per_m=0.0,connection_price_each=0.0,
                     insulation_price_per_m2=None,max_dp_Pa=None,solver_options=None):
    if fluid.lower()!='water':
        raise ValueError('Current validated fluid model: water only.')
    if heat_required<=0 or solar_irradiance<=0 or mdot_total<=0:
        raise ValueError('Heat, irradiance and mass flow must be positive.')
    if beta<0 or beta>90:
        raise ValueError('Collector tilt angle beta must be between 0 and 90 degrees.')
    if cover2_id and cover2_id==cover1_id:
        # Allowed: two identical covers are a legitimate design.
        pass
    if cover_gap_m<0:
        raise ValueError('Cover gap cannot be negative.')
    if connection_type not in ('below_plate','above_plate','in_line'):
        raise ValueError('Unknown connection type.')
    if solver_options is None:
        solver_options={'nx':16,'ny_per_gap':2,'max_outer':15}

    absorber=get_by_id(catalog('absorbers'),absorber_id)
    tube=get_by_id(catalog('tubes'),tube_id)
    back_ins_id = back_insulation_id if insulation_id is None else insulation_id
    back_insulation=get_by_id(catalog('insulation'),back_ins_id)
    side_insulation=get_by_id(catalog('insulation'),side_insulation_id)
    adhesive=get_by_id(catalog('adhesives'),adhesive_id)
    connections=get_by_id(catalog('connection_types'),connection_type)
    covers=cover_set(cover1_id,cover2_id)
    validate_tube(tube)

    do=_f(tube,'od_m'); di=_f(tube,'id_m')
    plate_k=_f(absorber,'thermal_conductivity_W_mK'); plate_delta=_f(absorber,'thickness_m')
    alpha=_f(absorber,'absorptivity'); eps_p=_f(absorber,'emissivity')
    tau,n_cover,cover_eps=effective_cover_properties(covers)
    M=len(covers)
    # The model applies tau**M, so pass the per-cover (geometric mean) value;
    # passing the product would apply the second cover twice.
    tau_per_cover=tau**(1.0/M)
    back_ins_k=_f(back_insulation,'thermal_conductivity_W_mK')
    back_ins_delta=_f(back_insulation,'default_thickness_m') if back_insulation_thickness_m is None else float(back_insulation_thickness_m)
    side_ins_k=_f(side_insulation,'thermal_conductivity_W_mK')
    side_ins_delta=_f(side_insulation,'default_thickness_m') if side_insulation_thickness_m is None else float(side_insulation_thickness_m)
    if back_ins_delta<=0 or side_ins_delta<=0: raise ValueError('Back and side insulation thicknesses must be positive.')
    adh_k=_f(adhesive,'thermal_conductivity_W_mK') if adhesive_k_W_mK is None else float(adhesive_k_W_mK)
    adh_delta=_f(adhesive,'default_thickness_m') if adhesive_thickness_m is None else float(adhesive_thickness_m)
    if adh_k<=0: raise ValueError('Adhesive thermal conductivity must be positive.')
    if adh_delta<=0: raise ValueError('Adhesive thickness must be positive for the current 2-D engineering model.')

    # Legacy BASIC only asked for ka. The thickness is retained here because
    # the modern 2-D model represents an explicit adhesive thermal resistance.
    # For in-line tubes the coupling is approximated as direct contact.
    model_adh_k=adh_k; model_adh_delta=adh_delta
    connection_note=connections['model_note']
    if connection_type=='in_line':
        model_adh_k=1.0e6
        model_adh_delta=1.0e-9

    if tube_price_per_m is None: tube_price_per_m=_f(tube,'cost_per_m')
    if absorber_price_per_m2 is None: absorber_price_per_m2=_f(absorber,'cost_per_m2')
    if cover_price_per_m2 is None: cover_price_per_m2=sum(_f(c,'cost_per_m2') for c in covers)
    if insulation_price_per_m2 is None: insulation_price_per_m2=_f(back_insulation,'cost_per_m2')

    feasible=[]; attempted=0
    for N in N_values:
        for w in spacing_values:
            attempted+=1
            try:
                r=solve_L1_for_q_fast(
                    mdot_total=mdot_total,N=N,w=w,q_target=heat_required,
                    di=di,do=do,plate_k=plate_k,plate_delta=plate_delta,
                    adhesive_k=model_adh_k,adhesive_delta=model_adh_delta,
                    insulation_k=back_ins_k,deltab=back_ins_delta,side_insulation_k=side_ins_k,side_insulation_thickness_m=side_ins_delta,
                    solar_I=solar_irradiance,Ta=ambient_temperature,Tfi=inlet_temperature,
                    wind=wind_speed,beta=beta,alpha=alpha,tau=tau_per_cover,n_cover=n_cover,eps_p=eps_p,M=M,cover_emissivity=cover_eps,
                    **solver_options)
            except Exception:
                continue
            if r['q_W']+0.2<heat_required: continue
            if r['L1_m']>max_L1 or r['L2_m']>max_L2: continue
            if r['flow_nonuniformity']>max_nonuniformity_percent/100.0: continue
            if max_dp_Pa is not None and r['dp_system_Pa']>max_dp_Pa: continue
            tube_length=N*r['L1_m']
            cost=(tube_length*tube_price_per_m + r['Ap_m2']*absorber_price_per_m2 +
                  r['Ap_m2']*cover_price_per_m2 + r['Ap_m2']*insulation_price_per_m2 +
                  2*r['L2_m']*header_price_per_m + N*connection_price_each)
            x=dict(r)
            x.update({'beta_deg':beta,'solar_irradiance_W_m2':solar_irradiance,'absorber_id':absorber_id,'tube_id':tube_id,'cover1_id':cover1_id,
                      'cover2_id':cover2_id or 'none','cover_count':M,'cover_gap_m':cover_gap_m,
                      'back_insulation_id':back_ins_id,'back_insulation_thickness_m':back_ins_delta,
                      'side_insulation_id':side_insulation_id,'side_insulation_thickness_m':side_ins_delta,
                      'adhesive_id':adhesive_id,'adhesive_k_W_mK':adh_k,'adhesive_thickness_m':adh_delta,
                      'connection_type':connection_type,'connection_name':connections['name'],
                      'tube_wall_m':0.5*(do-di),'tube_length_m':tube_length,'cost':cost})
            feasible.append(x)
    if not feasible:
        raise ValueError('No feasible design found. Relax constraints or enlarge the design space.')
    if objective=='min_cost': selected=min(feasible,key=lambda x:x['cost'])
    elif objective=='max_efficiency': selected=max(feasible,key=lambda x:x['efficiency'])
    elif objective=='min_area': selected=min(feasible,key=lambda x:x['Ap_m2'])
    else: raise ValueError('objective must be min_cost, max_efficiency or min_area')
    return {'selected':selected,'feasible_candidates':feasible,'attempted_candidates':attempted,
            'catalogs':{'absorber':absorber,'tube':tube,'cover1':covers[0],
                        'cover2':covers[1] if len(covers)>1 else None,'back_insulation':back_insulation,'side_insulation':side_insulation,
                        'adhesive':adhesive,'connection':connections},
            'model_notes':{'cover_model':'legacy-equivalent for one cover; engineering effective-property approximation for two covers',
                           'cover_gap_used':False,
                           'side_insulation_model':'Engineering extension: independent side k/thickness; original BASIC did not expose separate side insulation.',
                           'connection_model':connection_note}}


def pareto_candidates(feasible):
    front=[]
    for a in feasible:
        dominated=False
        for b in feasible:
            if b is a: continue
            if b['cost']<=a['cost'] and b['efficiency']>=a['efficiency'] and (b['cost']<a['cost'] or b['efficiency']>a['efficiency']):
                dominated=True; break
        if not dominated: front.append(a)
    return sorted(front,key=lambda x:x['cost'])


def format_report(result):
    r=result['selected']; c=result['catalogs']; lines=[]; add=lines.append
    add('COLLECTOR DESIGN ENGINE v1.5'); add('='*68); add('SELECTED DESIGN'); add('-'*68)
    add('Connection         : %s'%r['connection_name'])
    add('Tube               : %s'%c['tube'].get('name',''))
    add('Tube OD / ID       : %.3f / %.3f mm'%(1000*_f(c['tube'],'od_m'),1000*_f(c['tube'],'id_m')))
    add('Tube wall          : %.3f mm'%(1000*r['tube_wall_m']))
    add('Absorber           : %s'%c['absorber'].get('name',''))
    add('Back insulation    : %s, %.1f mm, k=%.4f W/m.K'%(c['back_insulation'].get('name',''),1000*r['back_insulation_thickness_m'],_f(c['back_insulation'],'thermal_conductivity_W_mK')))
    add('Side insulation    : %s, %.1f mm, k=%.4f W/m.K'%(c['side_insulation'].get('name',''),1000*r['side_insulation_thickness_m'],_f(c['side_insulation'],'thermal_conductivity_W_mK')))
    add('Adhesive           : %s, k=%.4f W/m.K'%(c['adhesive'].get('name',''),r['adhesive_k_W_mK']))
    add('Covers             : %d'%r['cover_count'])
    add('Cover 1            : %s'%c['cover1'].get('name',''))
    if c['cover2'] is not None: add('Cover 2            : %s'%c['cover2'].get('name',''))
    add('Cover gap          : %.1f mm'%(1000*r['cover_gap_m']))
    add(''); add('OPERATING / ENVIRONMENT'); add('Solar irradiance  : %.2f W/m2'%r['solar_irradiance_W_m2']); add('Collector tilt beta: %.2f deg from horizontal'%r['beta_deg']); add(''); add('GEOMETRY'); add('N tubes            : %d'%r['N']); add('Tube spacing       : %.3f m'%r['w_m']); add('Collector length   : %.3f m'%r['L1_m']); add('Collector width    : %.3f m'%r['L2_m']); add('Absorber area      : %.3f m2'%r['Ap_m2']); add('Total tube length  : %.3f m'%r['tube_length_m'])
    add(''); add('THERMAL'); add('Required heat      : %.2f W'%r['q_W']); add('Efficiency         : %.3f %%'%(100*r['efficiency'])); add('T plate max        : %.2f C'%r['Tmax_C']); add('T outlet mean      : %.2f C'%r['Tfo_mean_C']); add('Overall UL         : %.3f W/m2.K'%r['UL_W_m2K'])
    add(''); add('HYDRAULIC'); add('Pressure drop      : %.2f Pa'%r['dp_system_Pa']); add('Flow nonuniformity : %.3f %%'%(100*r['flow_nonuniformity'])); add('Pump power         : %.6f W'%r['pump_power_W'])
    add(''); add('ECONOMIC'); add('Entered/model cost  : %.4f'%r['cost'])
    add(''); add('STATUS / MODEL NOTES'); add('Thermal requirement: PASS'); add('Geometry limits    : PASS'); add('Flow uniformity    : PASS')
    if result['model_notes']['cover_gap_used'] is False: add('Cover gap          : recorded, not yet active in legacy top-loss equation')
    add('Cover model        : %s'%result['model_notes']['cover_model'])
    add('Side insulation    : independent engineering extension active')
    add('Feasible candidates: %d'%len(result['feasible_candidates']))
    return '\n'.join(lines)


if __name__=='__main__':
    res=design_collector(2090,20,15,760,2,mdot_total=.10,solver_options={'nx':12,'ny_per_gap':2,'max_outer':12})
    print(format_report(res))
