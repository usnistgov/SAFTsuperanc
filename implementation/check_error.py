import pandas
import glob
import json
import numpy as np
import PCSAFTsuperanc
import matplotlib.pyplot as plt
import teqp 

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
        
        def add_p(row):
            [tilderhoL, tilderhoV] = PCSAFTsuperanc.PCSAFTsuperanc_rhoLV(row['Ttilde'], m)
            sigma_m = 3e-10
            ek = 150
            T = row['Ttilde']*ek
            # N_A = 8.31446261815324/1.380649e-23
            # N_A = 6.022e23
            N_A = PCSAFTsuperanc.N_A
            rhoL, rhoV = [tilderho/(N_A*sigma_m**3) for tilderho in [tilderhoL, tilderhoV]]
            c = teqp.SAFTCoeffs()
            c.sigma_Angstrom = sigma_m*1e10
            c.epsilon_over_k = ek 
            c.m = m
            model = teqp.PCSAFTEOS([c])
            z = np.array([1.0])
            pL = rhoL*model.get_R(z)*T*(1+model.get_Ar01(T, rhoL, z))
            pV = rhoV*model.get_R(z)*T*(1+model.get_Ar01(T, rhoV, z))
            print(pL, pV, 'superanc densities')

            print([rhoL, rhoV], teqp.pure_VLE_T(model, T, rhoL, rhoV, 10))
            rhoL, rhoV = teqp.pure_VLE_T(model, T, rhoL, rhoV, 10)
            pL = rhoL*model.get_R(z)*T*(1+model.get_Ar01(T, rhoL, z))
            pV = rhoV*model.get_R(z)*T*(1+model.get_Ar01(T, rhoV, z))
            print(pL, pV, 'VLE densities')

            return pL, pV
            
        gp[['pL / Pa', 'pV / Pa']] = gp.apply(add_p, result_type='expand', axis=1)

        Ttilde_crit, Ttilde_min = PCSAFTsuperanc.get_Ttilde_crit_min(m)
        gp = gp[gp['Ttilde'] > Ttilde_min]

        # errL = np.abs(gp['rhotildeL(anc)']/gp['rhotildeL']-1)
        # print(m, 'L', np.max(errL))
        # errV = np.abs(gp['rhotildeV(anc)']/gp['rhotildeV']-1)
        # print(m, 'V', np.max(errV))

        # plt.plot(gp['rhotildeL(anc)'], gp['Ttilde'])
        # plt.plot(gp['rhotildeL'], gp['Ttilde'])
        # plt.show()

if __name__ == '__main__':
    pkl = 'all_VLE.pkl.xz'
    # rebuild_pickle(pkl)
    check_deviation(pkl)