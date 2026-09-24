# Tkinter GUI for Collector Design Engine v1.3
# Python 3.7.2 / standard library only.
import tkinter as tk
from tkinter import ttk, messagebox
from collector_design_engine_v1_2 import design_collector, catalog, format_report, pareto_candidates

class App(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self); self.title('Solar Collector Design Engine v1.3'); self.geometry('980x900'); self._build()
    def entry(self,parent,label,default,row,col=0,width=15):
        ttk.Label(parent,text=label).grid(row=row,column=col,sticky='w',padx=5,pady=3)
        e=ttk.Entry(parent,width=width); e.insert(0,str(default)); e.grid(row=row,column=col+1,sticky='w',padx=5,pady=3); return e
    def combo(self,parent,label,rows,row,col=0,width=28):
        ttk.Label(parent,text=label).grid(row=row,column=col,sticky='w',padx=5,pady=3)
        cb=ttk.Combobox(parent,values=[x['id'] for x in rows],state='readonly',width=width); cb.set(rows[0]['id']); cb.grid(row=row,column=col+1,sticky='w',padx=5,pady=3); return cb
    def _build(self):
        f=ttk.Frame(self,padding=8); f.pack(fill='both',expand=True)
        inp=ttk.LabelFrame(f,text='1. Design / environment'); inp.pack(fill='x')
        self.q=self.entry(inp,'Required heat Q [W]',2090,0,0); self.tfi=self.entry(inp,'Inlet temperature [C]',20,0,2); self.ta=self.entry(inp,'Ambient temperature [C]',15,1,0); self.I=self.entry(inp,'Solar irradiance [W/m2]',760,1,2); self.wind=self.entry(inp,'Wind speed [m/s]',2,2,0); self.beta=self.entry(inp,'Collector tilt beta [deg]',30,2,2); self.mdot=self.entry(inp,'Mass flow [kg/s]',0.10,3,0)
        mat=ttk.LabelFrame(f,text='2. Materials'); mat.pack(fill='x',pady=6)
        self.abs_rows=catalog('absorbers'); self.tube_rows=catalog('tubes'); self.cover_rows=catalog('covers'); self.ins_rows=catalog('insulation'); self.adh_rows=catalog('adhesives'); self.conn_rows=catalog('connection_types')
        self.abs=self.combo(mat,'Absorber',self.abs_rows,0); self.tube=self.combo(mat,'Tube',self.tube_rows,1)
        self.ins=self.combo(mat,'Back insulation',self.ins_rows,2); self.adh=self.combo(mat,'Adhesive',self.adh_rows,3)
        self.conn=self.combo(mat,'Tube/plate connection',self.conn_rows,4)
        ttk.Label(mat,text='  1 below plate / 2 above plate / 3 in-line').grid(row=4,column=2,sticky='w')
        cov=ttk.LabelFrame(f,text='3. Covers'); cov.pack(fill='x',pady=6)
        self.cover_count=tk.IntVar(value=1); ttk.Label(cov,text='Number of covers').grid(row=0,column=0,sticky='w',padx=5,pady=3); ttk.Combobox(cov,textvariable=self.cover_count,values=[1,2],state='readonly',width=8).grid(row=0,column=1,sticky='w')
        self.c1=self.combo(cov,'Cover 1',self.cover_rows,1); self.c2=self.combo(cov,'Cover 2',self.cover_rows,2); self.cover_gap=self.entry(cov,'Gap between covers [mm]',12,3); self.c2.configure(state='disabled'); self.cover_count.trace('w',self._toggle_cover2)
        eng=ttk.LabelFrame(f,text='4. Engineering dimensions'); eng.pack(fill='x',pady=6)
        self.ins_th=self.entry(eng,'Insulation thickness [mm]',40,0,0); self.adh_th=self.entry(eng,'Adhesive thickness [mm]',0.5,0,2); self.adh_k=self.entry(eng,'Adhesive k [W/m.K]',0.2,1,0); self.maxL1=self.entry(eng,'Max collector length [m]',3.5,1,2); self.maxL2=self.entry(eng,'Max collector width [m]',1.5,2,0); self.maxnu=self.entry(eng,'Max flow nonuniformity [%]',5,2,2)
        cost=ttk.LabelFrame(f,text='5. Cost / constraints'); cost.pack(fill='x',pady=6)
        self.pt=self.entry(cost,'Tube price [currency/m]',1,0,0); self.pa=self.entry(cost,'Absorber price [currency/m2]',1,0,2); self.pg=self.entry(cost,'Cover price [currency/m2]',1,1,0); self.pi=self.entry(cost,'Insulation price [currency/m2]',1,1,2); self.ph=self.entry(cost,'Header price [currency/m]',0,2,0); self.pj=self.entry(cost,'Connection price each',0,2,2); self.maxdp=self.entry(cost,'Max pressure drop [Pa]',5000,3,0)
        opt=ttk.Frame(f); opt.pack(fill='x',pady=6); ttk.Label(opt,text='Objective').pack(side='left'); self.obj=ttk.Combobox(opt,values=['min_cost','max_efficiency','min_area'],state='readonly',width=18); self.obj.set('min_cost'); self.obj.pack(side='left',padx=8); ttk.Button(opt,text='OPTIMIZE',command=self.run).pack(side='left',padx=15)
        out=ttk.LabelFrame(f,text='6. Engineering result'); out.pack(fill='both',expand=True); self.text=tk.Text(out,font=('Consolas',9),wrap='none'); self.text.pack(fill='both',expand=True)
    def _toggle_cover2(self,*args):
        if self.cover_count.get()==2: self.c2.configure(state='readonly')
        else: self.c2.configure(state='disabled')
    def run(self):
        try:
            c2=self.c2.get() if self.cover_count.get()==2 else None
            r=design_collector(float(self.q.get()),float(self.tfi.get()),float(self.ta.get()),float(self.I.get()),float(self.wind.get()),float(self.beta.get()),
                absorber_id=self.abs.get(),tube_id=self.tube.get(),cover1_id=self.c1.get(),cover2_id=c2,cover_gap_m=float(self.cover_gap.get())/1000.0,
                insulation_id=self.ins.get(),insulation_thickness_m=float(self.ins_th.get())/1000.0,adhesive_id=self.adh.get(),adhesive_k_W_mK=float(self.adh_k.get()),adhesive_thickness_m=float(self.adh_th.get())/1000.0,
                connection_type=self.conn.get(),mdot_total=float(self.mdot.get()),objective=self.obj.get(),
                tube_price_per_m=float(self.pt.get()),absorber_price_per_m2=float(self.pa.get()),cover_price_per_m2=float(self.pg.get()),insulation_price_per_m2=float(self.pi.get()),header_price_per_m=float(self.ph.get()),connection_price_each=float(self.pj.get()),max_dp_Pa=float(self.maxdp.get()),
                max_L1=float(self.maxL1.get()),max_L2=float(self.maxL2.get()),max_nonuniformity_percent=float(self.maxnu.get()),N_values=range(8,21),spacing_values=(.06,.08,.10,.12,.14,.16),solver_options={'nx':12,'ny_per_gap':2,'max_outer':12})
            self.text.delete('1.0','end'); self.text.insert('1.0',format_report(r)); self.text.insert('end','\n\nPARETO CANDIDATES\n------------------\n')
            for x in pareto_candidates(r['feasible_candidates']): self.text.insert('end','N=%d  w=%.3f  cost=%.3f  eta=%.3f%%  conn=%s\n'%(x['N'],x['w_m'],x['cost'],100*x['efficiency'],x['connection_type']))
        except Exception as e: messagebox.showerror('Design error',str(e))

if __name__=='__main__': App().mainloop()
