# Improved Tkinter GUI for Solar Collector Design Engine v1.5
# Python 3.7.2 / standard library only.
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
from collector_design_engine_v1_3 import design_collector, catalog, format_report, pareto_candidates

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Solar Collector Design Engine v1.5')
        self.geometry('1120x820'); self.minsize(980,700)
        self._make_style(); self._build()
    def _make_style(self):
        s=ttk.Style(self)
        try: s.theme_use('clam')
        except Exception: pass
        s.configure('Title.TLabel',font=('Segoe UI',15,'bold'))
        s.configure('Hint.TLabel',foreground='#555555')
        s.configure('Run.TButton',font=('Segoe UI',10,'bold'))
    def field(self,p,label,default,row,col,width=14):
        ttk.Label(p,text=label).grid(row=row,column=col,sticky='w',padx=6,pady=4)
        e=ttk.Entry(p,width=width); e.insert(0,str(default)); e.grid(row=row,column=col+1,sticky='ew',padx=6,pady=4); return e
    def combo(self,p,label,rows,row,col,width=25):
        ttk.Label(p,text=label).grid(row=row,column=col,sticky='w',padx=6,pady=4)
        c=ttk.Combobox(p,values=[x['id'] for x in rows],state='readonly',width=width)
        if rows: c.current(0)
        c.grid(row=row,column=col+1,sticky='ew',padx=6,pady=4); return c
    def _build(self):
        root=ttk.Frame(self,padding=10); root.pack(fill='both',expand=True)
        ttk.Label(root,text='Solar Collector Design Engine',style='Title.TLabel').pack(anchor='w')
        ttk.Label(root,text='Legacy BASIC correlations preserved; engineering extensions are reported separately.',style='Hint.TLabel').pack(anchor='w',pady=(0,8))
        nb=ttk.Notebook(root); nb.pack(fill='both',expand=True)
        self.tab1=ttk.Frame(nb,padding=8); self.tab2=ttk.Frame(nb,padding=8); self.tab3=ttk.Frame(nb,padding=8)
        nb.add(self.tab1,text='Conditions'); nb.add(self.tab2,text='Materials and geometry'); nb.add(self.tab3,text='Cost and constraints')
        self._conditions(); self._materials(); self._costs()
        bottom=ttk.Frame(root); bottom.pack(fill='x',pady=(8,0))
        self.status=ttk.Label(bottom,text='Ready.',style='Hint.TLabel'); self.status.pack(side='left')
        ttk.Button(bottom,text='Open data folder',command=self.open_data).pack(side='right',padx=4)
        ttk.Button(bottom,text='Save report',command=self.save_report).pack(side='right',padx=4)
        ttk.Button(bottom,text='OPTIMIZE',style='Run.TButton',command=self.run).pack(side='right',padx=4)
        out=ttk.LabelFrame(root,text='Results'); out.pack(fill='both',expand=True,pady=(8,0))
        self.text=tk.Text(out,font=('Consolas',9),wrap='none'); self.text.pack(side='left',fill='both',expand=True)
        sb=ttk.Scrollbar(out,orient='vertical',command=self.text.yview); sb.pack(side='right',fill='y'); self.text.configure(yscrollcommand=sb.set)
    def _conditions(self):
        p=ttk.LabelFrame(self.tab1,text='Design point'); p.pack(fill='x');
        self.q=self.field(p,'Required heat Q [W]',2090,0,0); self.tfi=self.field(p,'Inlet temperature [°C]',20,0,2)
        self.ta=self.field(p,'Ambient temperature [°C]',15,1,0); self.I=self.field(p,'Solar irradiance [W/m²]',760,1,2)
        self.wind=self.field(p,'Wind speed [m/s]',2,2,0); self.beta=self.field(p,'Tilt angle β [deg]',30,2,2)
        self.mdot=self.field(p,'Mass flow [kg/s]',0.10,3,0)
        ttk.Label(p,text='Irradiance is entered directly on the collector plane; no solar-position model is used.',style='Hint.TLabel').grid(row=4,column=0,columnspan=4,sticky='w',padx=6,pady=6)
    def _materials(self):
        self.rows={n:catalog(n) for n in ['absorbers','tubes','covers','insulation','adhesives','connection_types']}
        p=ttk.LabelFrame(self.tab2,text='Material selection'); p.pack(fill='x')
        self.abs=self.combo(p,'Absorber',self.rows['absorbers'],0,0); self.tube=self.combo(p,'Tube',self.rows['tubes'],1,0)
        self.back_ins=self.combo(p,'Back insulation',self.rows['insulation'],2,0); self.side_ins=self.combo(p,'Side insulation',self.rows['insulation'],3,0); self.adh=self.combo(p,'Adhesive',self.rows['adhesives'],4,0)
        self.conn=self.combo(p,'Tube/plate connection',self.rows['connection_types'],5,0)
        self.back_info=ttk.Label(p,text='',style='Hint.TLabel'); self.back_info.grid(row=2,column=2,columnspan=2,sticky='w',padx=6)
        self.side_info=ttk.Label(p,text='',style='Hint.TLabel'); self.side_info.grid(row=3,column=2,columnspan=2,sticky='w',padx=6)
        self.back_ins.bind('<<ComboboxSelected>>',self.update_insulation_info); self.side_ins.bind('<<ComboboxSelected>>',self.update_insulation_info)
        c=ttk.LabelFrame(self.tab2,text='Covers and layers'); c.pack(fill='x',pady=8)
        self.cover_count=tk.IntVar(value=1); ttk.Label(c,text='Number of covers').grid(row=0,column=0,sticky='w',padx=6,pady=4)
        ttk.Combobox(c,textvariable=self.cover_count,values=[1,2],state='readonly',width=8).grid(row=0,column=1,sticky='w',padx=6,pady=4)
        self.c1=self.combo(c,'Cover 1',self.rows['covers'],1,0); self.c2=self.combo(c,'Cover 2',self.rows['covers'],2,0)
        self.cover_gap=self.field(c,'Gap between covers [mm]',12,3,0); self.c2.configure(state='disabled'); self.cover_count.trace('w',self.toggle2)
        self.back_ins_th=self.field(c,'Back insulation thickness [mm]',40,4,0); self.side_ins_th=self.field(c,'Side insulation thickness [mm]',40,5,0); self.adh_th=self.field(c,'Adhesive thickness [mm]',0.5,4,2); self.adh_k=self.field(c,'Adhesive k [W/m·K]',0.2,6,0)
        self.update_insulation_info()
    def _costs(self):
        p=ttk.LabelFrame(self.tab3,text='Prices'); p.pack(fill='x')
        self.pt=self.field(p,'Tube [currency/m]',1,0,0); self.pa=self.field(p,'Absorber [currency/m²]',1,0,2); self.pg=self.field(p,'Cover [currency/m²]',1,1,0); self.pi=self.field(p,'Insulation [currency/m²]',1,1,2); self.ph=self.field(p,'Header [currency/m]',0,2,0); self.pj=self.field(p,'Connection each',0,2,2)
        q=ttk.LabelFrame(self.tab3,text='Constraints and objective'); q.pack(fill='x',pady=8)
        self.maxL1=self.field(q,'Max length L1 [m]',3.5,0,0); self.maxL2=self.field(q,'Max width L2 [m]',1.5,0,2); self.maxnu=self.field(q,'Max nonuniformity [%]',5,1,0); self.maxdp=self.field(q,'Max pressure drop [Pa]',5000,1,2)
        ttk.Label(q,text='Objective').grid(row=2,column=0,sticky='w',padx=6,pady=4); self.obj=ttk.Combobox(q,values=['min_cost','max_efficiency','min_area'],state='readonly',width=18); self.obj.set('min_cost'); self.obj.grid(row=2,column=1,sticky='w',padx=6,pady=4)
    def update_insulation_info(self,*a):
        for cb,label in ((self.back_ins,self.back_info),(self.side_ins,self.side_info)):
            try:
                r=next(x for x in self.rows['insulation'] if x['id']==cb.get())
                label.config(text='k=%.4f W/m.K | default=%.1f mm | status=%s'%(float(r['thermal_conductivity_W_mK']),1000*float(r['default_thickness_m']),r.get('data_status','')))
            except Exception:
                label.config(text='')
    def toggle2(self,*a): self.c2.configure(state='readonly' if self.cover_count.get()==2 else 'disabled')
    def open_data(self):
        path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')
        try: os.startfile(path)
        except Exception: messagebox.showinfo('Data folder',path)
    def save_report(self):
        txt=self.text.get('1.0','end').strip()
        if not txt: return
        fn=filedialog.asksaveasfilename(defaultextension='.txt',filetypes=[('Text','*.txt')])
        if fn:
            with open(fn,'w',encoding='utf-8') as f: f.write(txt)
    def run(self):
        try:
            self.status.config(text='Calculating...'); self.update_idletasks()
            c2=self.c2.get() if self.cover_count.get()==2 else None
            r=design_collector(float(self.q.get()),float(self.tfi.get()),float(self.ta.get()),float(self.I.get()),float(self.wind.get()),float(self.beta.get()),
                absorber_id=self.abs.get(),tube_id=self.tube.get(),cover1_id=self.c1.get(),cover2_id=c2,cover_gap_m=float(self.cover_gap.get())/1000,
                back_insulation_id=self.back_ins.get(),back_insulation_thickness_m=float(self.back_ins_th.get())/1000,side_insulation_id=self.side_ins.get(),side_insulation_thickness_m=float(self.side_ins_th.get())/1000,adhesive_id=self.adh.get(),adhesive_k_W_mK=float(self.adh_k.get()),adhesive_thickness_m=float(self.adh_th.get())/1000,
                connection_type=self.conn.get(),mdot_total=float(self.mdot.get()),objective=self.obj.get(),tube_price_per_m=float(self.pt.get()),absorber_price_per_m2=float(self.pa.get()),cover_price_per_m2=float(self.pg.get()),insulation_price_per_m2=float(self.pi.get()),header_price_per_m=float(self.ph.get()),connection_price_each=float(self.pj.get()),max_dp_Pa=float(self.maxdp.get()),max_L1=float(self.maxL1.get()),max_L2=float(self.maxL2.get()),max_nonuniformity_percent=float(self.maxnu.get()),N_values=range(8,21),spacing_values=(.06,.08,.10,.12,.14,.16),solver_options={'nx':12,'ny_per_gap':2,'max_outer':12})
            self.text.delete('1.0','end'); self.text.insert('1.0',format_report(r)); self.text.insert('end','\n\nPARETO CANDIDATES\n------------------\n')
            for x in pareto_candidates(r.get('feasible_candidates',[])): self.text.insert('end','N=%d  w=%.3f  cost=%.3f  eta=%.3f%%  conn=%s\n'%(x['N'],x['w_m'],x['cost'],100*x['efficiency'],x['connection_type']))
            self.status.config(text='Calculation completed.')
        except Exception as e:
            self.status.config(text='Error.'); messagebox.showerror('Design error',str(e))
if __name__=='__main__': App().mainloop()
