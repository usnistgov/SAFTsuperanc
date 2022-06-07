import ChebTools, numpy as np, json

class Evaluator:
    mmin, mmax = 1, 64
    ymin, ymax = 1/mmax, 1/mmin # y=1/m
    expsL, expsV = [], []
    N = 64 # degree in 1/m direction
    k = np.arange(0,N+1,1)
    ynodes = (ymax-ymin)/2*np.cos(k*np.pi/N) + (ymin+ymax)/2
    mnodes = 1/ynodes
    print(mnodes)

    def __init__(self):
        print(self.mnodes)
        """ Unpack the expansions """
        for m in self.mnodes:
            L, V = [], []
            for ex in json.load(open(f'bld/PCSAFT_VLE_m{m:0.6f}_expansions.json')):
                L.append(ChebTools.ChebyshevExpansion(ex['coefL'],ex['xmin'],ex['xmax']))
                V.append(ChebTools.ChebyshevExpansion(ex['coefV'],ex['xmin'],ex['xmax']))
            self.expsL.append(ChebTools.ChebyshevCollection(L))
            self.expsV.append(ChebTools.ChebyshevCollection(V))

    def get_node_vals(self, *, Theta):
        rhoLfvals = np.array([ex(Theta) for ex in self.expsL])
        rhoVfvals = np.array([ex(Theta) for ex in self.expsV])
        return rhoLfvals, rhoVfvals
        
    def __call__(self, *, Theta, m):
        """ Call the function to get densities """
        rhoLfvals, rhoVfvals = self.get_node_vals(Theta=Theta)
        ce_rhoL = ChebTools.factoryfDCT(self.N, rhoLfvals, self.ymin, self.ymax)
        tilderhoL = ce_rhoL.y(1/m)
        ce_rhoV = ChebTools.factoryfDCT(self.N, rhoVfvals, self.ymin, self.ymax)
        tilderhoV = ce_rhoV.y(1/m)
        return (tilderhoL, tilderhoV)

e = Evaluator()

m = 4.0
print(e(Theta=0.95, m=m))

import matplotlib.pyplot as plt
for Theta in [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]:
    plt.plot(1/e.mnodes, e.get_node_vals(Theta=Theta)[1], 'o-', label=Theta)
plt.legend()
plt.show()

import pandas
j = json.load(open(f'bld/PCSAFT_VLE_m{m:06f}.json'))
arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
Ttildec = j['Ttildec']

# From the empirical fit
c = [ 0.37627892, -2.20078778 ]
Tredmin = np.exp(c[1])*np.power(m, c[0])
Ttildemin = Tredmin*j['Ttildec']

df = pandas.DataFrame({k:j[k] for k in arrays})
df['Theta'] = (df['Ttilde']-Ttildemin)/(Ttildec-Ttildemin)
df.to_csv(f'{m:0.6f}.csv', index=False)