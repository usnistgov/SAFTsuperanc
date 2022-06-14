import pandas
import glob
import json
import numpy as np
import PCSAFTsuperanc
import matplotlib.pyplot as plt

def rebuild_pickle(pkl):
    dfs = []
    for fname in glob.glob('../output/PCSAFT_VLE_m*.json'):
        if 'expansions' in fname:
            continue
        print(fname)
        j = json.load(open(fname))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        df['m'] = j['m']
        dfs.append(df)
    df = pandas.concat(dfs,sort=False).sort_values(by=['m','Ttilde'])
    df.to_pickle(pkl)

def check_deviation(pkl):
    df = pandas.read_pickle(pkl).sort_values(by=['m', 'Ttilde'])

    for m, gp in df.groupby('m'):
        gp = gp.copy()
        def add_anc(row):
            return PCSAFTsuperanc.PCSAFTsuperanc_rhoLV(row['Ttilde'], m)
        gp[['rhotildeL(anc)', 'rhotildeV(anc)']] = gp.apply(add_anc, result_type='expand', axis=1)
        Ttilde_crit, Ttilde_min = PCSAFTsuperanc.get_Ttilde_crit_min(m)
        gp = gp[gp['Ttilde'] > Ttilde_min]

        errL = np.abs(gp['rhotildeL(anc)']/gp['rhotildeL']-1)
        print(m, 'L', np.max(errL))
        errV = np.abs(gp['rhotildeV(anc)']/gp['rhotildeV']-1)
        print(m, 'V', np.max(errV))

        # plt.plot(gp['rhotildeL(anc)'], gp['Ttilde'])
        # plt.plot(gp['rhotildeL'], gp['Ttilde'])
        # plt.show()

if __name__ == '__main__':
    pkl = 'all_VLE.pkl.xz'
    rebuild_pickle(pkl)
    check_deviation(pkl)