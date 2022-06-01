import glob
import pandas, json, matplotlib.pyplot as plt, numpy as np
import scipy.interpolate
plt.style.use('mystyle.mplstyle')

root = 'bld'

def plot_critical_curve():
    j = pandas.DataFrame(json.load(open(f'{root}/PCSAFT_crit_pts_interpolation.json')))
    plt.plot(j['1/m'], j['Ttilde'])
    plt.gca().set(xlabel='1/m', ylabel=r'$\widetilde{T}_{\rm crit}$')
    plt.show()

def get_fnames():
    fnames = glob.glob(f'{root}/PCSAFT_VLE_m*.json')
    ms, names = zip(*sorted([(json.load(open(fname))['m'], fname) for fname in fnames if '_expansions' not in fname]))
    return names

def plot_all_VLE():

    fig, ax = plt.subplots(1,1,figsize=(3.3,3))
    for f in get_fnames():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        if m > 64 and m != 80:
            continue
        Ttildemin = np.min(df['Ttilde'])
        
        Tred = df['Ttilde']/j['Ttildec']
        line, = plt.plot(df['rhotildeL']/j['rhotildec'], Tred, label=f'$m$: {m}')
        plt.plot(df['rhotildeV']/j['rhotildec'], Tred, color=line.get_color())
    plt.axvline(1.0, dashes=[1,1])
    plt.gca().set(
        xlabel=r'$\widetilde{\rho}/\widetilde{\rho}_{\rm crit}$', 
        ylabel=r'$\widetilde{T}/\widetilde{T}_{\rm crit}$'
    )
    plt.legend(loc='best', ncol=3, fontsize=4)
    plt.tight_layout(pad=0.2)
    plt.savefig('PCSAFT_all_VLE.pdf')
    plt.close()

def plot_Tmin():
    ooo = []
    fig, ax = plt.subplots(1,1,figsize=(3.3,3))
    for f in get_fnames():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        Tred = df['Ttilde']/j['Ttildec']
        line, = plt.plot(df['rhotildeL']/j['rhotildeV'], Tred, label=f'$m$: {m}')
        try:
            Tmin_red = float(scipy.interpolate.interp1d(df['rhotildeL']/df['rhotildeV'], Tred)(1e20))
            ooo.append({'m': j['m'], 'Tmin_red': Tmin_red})
        except:
            pass
    plt.legend(loc='best')
    plt.xscale('log')
    plt.xlim(1, 1e40)
    plt.gca().set(
        xlabel=r"$\widetilde{\rho}'/\widetilde{\rho}''$",
        ylabel=r'$\widetilde{T}/\widetilde{T}_{\rm crit}$'
    )
    plt.savefig('PCSAFT_Tmin_VLE.pdf')
    plt.close()

    fig, ax = plt.subplots(1,1,figsize=(3.3,2))
    df = pandas.DataFrame(ooo)
    plt.plot(df['m'], df['Tmin_red'], 'o')
    x = np.array(df['m']); y = np.array(df['Tmin_red'])
    c = np.polyfit(np.log(x), np.log(y), 1)
    m = np.geomspace(1, 100)
    plt.plot(m, np.exp(np.polyval(c, np.log(m))), dashes=[2, 2])
    plt.plot(m, np.exp(c[1])*m**c[0], dashes=[3, 1, 1, 1])
    # Fit of the form
    # ln(Ttildemin) = c0*ln(m) + c1
    # ...
    # Ttildemin = exp(ln(m**(c0)) + c1) = exp(c1)*m**(c0)
    
    # plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$m$')
    plt.ylabel(r'$\widetilde{T}_{{\rm red}, \min}$')
    plt.tight_layout(pad=0.2)
    plt.savefig('Tmin_m.pdf')
    plt.close()
    
def plot_normalized_VLE():
    fig, ax = plt.subplots(1,1,figsize=(3.3,3))
    for f in get_fnames():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        if m > 80:
            continue
        
        # Ttildemin = np.min(df['Ttilde']) # from the data

        # From the empirical fit
        c = [ 0.37627892, -2.20078778 ]
        Tredmin = np.exp(c[1])*np.power(m, c[0])
        Ttildemin = Tredmin*j['Ttildec']

        Theta = (df['Ttilde']-Ttildemin)/(j['Ttildec']-Ttildemin)
        
        line, = plt.plot(df['rhotildeL']/j['rhotildec'], Theta, label=f'$m$: {m}')
        plt.plot(df['rhotildeV']/j['rhotildec'], Theta, color=line.get_color())
    plt.axvline(1.0, dashes=[1,1])
    plt.gca().set(
        xlabel=r'$\widetilde{\rho}/\widetilde{\rho}_{\rm crit}$', 
        ylabel=r'$\Theta=(\widetilde{T}-\widetilde{T}_{\rm min})/(\widetilde{T}_{\rm crit}-\widetilde{T}_{\rm min})$'
    )
    plt.ylim(0, 1)
    plt.legend(loc='best', fontsize=4, ncol=3)
    plt.tight_layout(pad=0.2)
    plt.savefig('PCSAFT_normalized_VLE.pdf')
    plt.close()

plot_all_VLE()
plot_normalized_VLE()
plot_Tmin()
# plot_critical_curve()