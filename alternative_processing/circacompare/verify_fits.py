from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import t as tdist, f as fdist, false_discovery_control

base = Path("alternative_processing/circacompare/calls")
tables = Path("results/tables")
layout = pd.read_csv(base/"sample_hours.csv")
summary = pd.read_csv(base/"pair_summary.csv")
panel = {}
checks = []
for d, s in layout.groupby("dataset", sort=False):
    x = pd.read_csv(tables / (("wt_cartilage" if d=="cartilage_wt" else d) + "_circ.csv"))
    q = pd.read_csv(base/f"single_{d}.csv")
    assert q.ensembl_gene_id.equals(x["CycID" if d=="cartilage_wt" else "Row.names"])
    assert q.error.isna().all(), (d, q.loc[q.error.notna(), ["symbol", "error"]])
    y = x[s["sample"]].to_numpy()
    t = s.elapsed_or_ZT_hour.to_numpy()*2*np.pi/24
    design = np.column_stack((np.ones(len(t)), np.cos(t), np.sin(t)))
    coef = np.linalg.lstsq(design, y.T, rcond=None)[0].T
    amplitude = np.hypot(coef[:,1], coef[:,2])
    unit = coef[:,1:]/amplitude[:,None]
    df = len(t)-3
    rss = np.sum((y-coef@design.T)**2, axis=1)
    var = rss/df
    cov = np.linalg.inv(design.T @ design)
    amp_var = np.einsum('ij,jk,ik->i', unit, cov[1:,1:], unit)
    tangent = np.column_stack((-unit[:,1],unit[:,0]))/amplitude[:,None]
    phi_var = np.einsum('ij,jk,ik->i', tangent, cov[1:,1:], tangent)
    expected_p = 2*tdist.sf(amplitude/np.sqrt(amp_var*var),df)
    expected_q = false_discovery_control(expected_p)
    assert np.array_equal(expected_p < .05,q.rhythmic_p < .05)
    assert np.array_equal(expected_q < .05,q.rhythmic_q < .05)
    joint_p = fdist.sf(((y.var(axis=1)*len(t)-rss)/2)/var,2,df)
    np.testing.assert_allclose(q.harmonic_F_p,joint_p,atol=1e-11)
    phase = np.arctan2(coef[:,2],coef[:,1])*24/(2*np.pi)
    phase_err = np.max(np.abs((phase-q.peak_time_hours+12)%24-12))
    assert phase_err < .02
    panel[d] = pd.DataFrame(dict(mesor=coef[:,0],amplitude=amplitude,phase=phase,
                                 mesor_var=cov[0,0], amplitude_var=amp_var,phase_var=phi_var*(24/(2*np.pi))**2,
                                 rss=rss,df=df,p=q.rhythmic_p.to_numpy(),q=q.rhythmic_q.to_numpy()),index=q.ensembl_gene_id)
    checks.append((d,len(q),q.clock_shift_h.ne(0).sum(),max(abs(expected_p-q.rhythmic_p)),phase_err))
for first,second in combinations(panel,2):
    l = f"{first}_vs_{second}"
    pair = pd.read_csv(base/f"pair_{l}.csv").set_index("ensembl_gene_id")
    a,b = panel[first],panel[second]
    assert set(pair.index) == (set(a.index[a.p < .05]) & set(b.index[b.p < .05]))
    a,b = a.loc[pair.index],b.loc[pair.index]
    assert np.array_equal(pair.both_rhythmic_q05,(a.q<.05)&(b.q<.05))
    assert pair.error.isna().all(), (l, pair.loc[pair.error.notna(),"error"])
    df = a.df+b.df
    pooled = (a.rss+b.rss)/df
    errors = []
    for metric in ["mesor","amplitude","phase"]:
        diff = b[metric]-a[metric]
        if metric=="phase": diff=(diff+12)%24-12
        var=pooled*(a[metric+"_var"]+b[metric+"_var"])
        expected_p=2*tdist.sf(abs(diff)/np.sqrt(var),df)
        errors.append(max(abs(expected_p-pair[metric+"_p"])))
        error = pair["phase_delta_h" if metric=="phase" else metric+"_delta"]-diff
        assert (abs(error)/np.sqrt(var)).max()<.002, (l,metric,(abs(error)/np.sqrt(var)).max())
        for family,keep in [("nominal_set",np.ones(len(pair),dtype=bool)),("both_rhythmic_set",pair.both_rhythmic_q05)]:
            q=false_discovery_control(expected_p[keep])
            row=summary[(summary["first"]==first)&(summary["second"]==second)&(summary["family"]==family)].iloc[0]
            assert np.array_equal(q<.05,pair.loc[keep,metric+"_q_"+family]<.05),(l,metric,family)
            assert (q<.05).sum()==row[metric+"_q05"]
    checks.append((l,len(pair),pair.clock_shift_h.ne(0).sum(),max(errors),max(abs((b.phase-a.phase-pair.phase_delta_h+12)%24-12))))
print(pd.DataFrame(checks,columns=["fits","n","retries","max_abs_p_error","max_phase_error_hours"]).to_string(index=False))

