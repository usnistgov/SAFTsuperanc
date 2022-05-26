import glob
import pandas, json, matplotlib.pyplot as plt, numpy as np

def plot_critical_curve():
    j = pandas.DataFrame(json.load(open('bld/PCSAFT_crit_pts_interpolation.json')))
    plt.plot(j['1/m'], j['Ttilde'])
    plt.gca().set(xlabel='1/m', ylabel=r'$\widetilde{T}_{\rm crit}$')
    plt.show()

def plot_all_VLE():
    fnames = glob.glob('bld/PCSAFT_VLE_m*.json')
    ms, names = zip(*sorted([(json.load(open(fname))['m'], fname) for fname in fnames]))
    for f in names:
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        Ttildemin = np.min(df['Ttilde'])
        Theta = (df['Ttilde']-Ttildemin)/(j['Ttildec']-Ttildemin)
        m = j['m']
        print(m, Ttildemin, Ttildemin/j['Ttildec'])
        line, = plt.plot(df['rhotildeL']/j['rhotildec'], Theta, label=f'$m$: {m}')
        plt.plot(df['rhotildeV']/j['rhotildec'], Theta, color=line.get_color())
    plt.axvline(1.0, dashes=[1,1])
    plt.gca().set(
        xlabel=r'$\widetilde{\rho}/\widetilde{\rho}_{\rm crit}$', 
        ylabel=r'$\Theta=(\widetilde{T}-\widetilde{T}_{\rm min})/(\widetilde{T}_{\rm crit}-\widetilde{T}_{\rm min})$'
    )
    plt.legend(loc='best')
    plt.savefig('PCSAFT_normalized_VLE.pdf')
    plt.show()

plot_critical_curve()
plot_all_VLE()