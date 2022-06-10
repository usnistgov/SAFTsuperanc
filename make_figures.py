import glob
import pandas, json, matplotlib.pyplot as plt, numpy as np
import scipy.interpolate
import re
plt.style.use('classic')
plt.style.use('mystyle.mplstyle')

root = 'output'

def plot_critical_curve():
    j = pandas.DataFrame(json.load(open(f'{root}/PCSAFT_crit_pts_interpolation.json')))
    fig, axes = plt.subplots(2,1,figsize=(3.3, 4), sharex=True)
    axes[0].plot(j['1/m'], j['Ttilde'])
    axes[1].plot(j['1/m'], j['rhotilde'])
    axes[0].set(ylabel=r'$\widetilde{T}_{\rm crit}$')
    axes[1].set(xlabel='$1/m$', ylabel=r'$\widetilde{\rho}_{\rm crit}$')
    plt.xscale('log')
    axes[1].set_yscale('log')
    plt.tight_layout(pad=0.2)
    plt.savefig('critical_values.pdf')
    plt.show()

def plot_critical_curvedev():
    j = pandas.DataFrame(json.load(open(f'{root}/PCSAFT_crit_pts_check.json')))
    fig, axes = plt.subplots(2,1,figsize=(3.3, 4), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(j['1/m'], j['Ttilde_tab'])
    axes[1].plot(j['1/m'], (j['Ttilde_tab']/j['Ttilde_fit']-1))

    axes[0].set(ylabel=r'$\widetilde{T}_{\rm crit}$')
    axes[1].set(xlabel='$1/m$', ylabel=r'$\widetilde{T}_{\rm crit,tab}/\widetilde{T}_{\rm crit,fit}-1$')
    plt.xscale('log')
    plt.tight_layout(pad=0.2)
    plt.savefig('Tcritical_values_dev.pdf')
    plt.show()

    fig, axes = plt.subplots(2,1,figsize=(3.3, 4), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(j['1/m'], j['rhotilde_tab'])
    axes[1].plot(j['1/m'], (j['rhotilde_tab']/j['rhotilde_fit']-1))

    axes[0].set(ylabel=r'$\widetilde{\rho}_{\rm crit}$')
    axes[1].set(xlabel='$1/m$', ylabel=r'$\widetilde{\rho}_{\rm crit,tab}/\widetilde{\rho}_{\rm crit,fit}-1$')
    plt.xscale('log')
    axes[0].set_yscale('log')
    plt.tight_layout(pad=0.2)
    plt.savefig('rhocritical_values_dev.pdf')
    plt.show()

def get_fnames():
    fnames = glob.glob(f'{root}/PCSAFT_VLE_m*.json')
    ms, names = zip(*sorted([(json.load(open(fname))['m'], fname) for fname in fnames if '_expansions' not in fname]))
    return names

def get_fnames_edges():
    Wedge_files = glob.glob(f'{root}/Wedges_pass*.json')
    def get_pass(f):
        match = re.search(r'Wedges_pass([0-9]+).json', f).group(1)
        return int(match)
    Wedge_filess = sorted(Wedge_files, key=get_pass)
    fname = Wedge_filess[-1]
    import os
    fnames = [f'{root}/PCSAFT_VLE_m{1/f:0.6f}.json' for f in json.load(open(fname))["Wedges"]]
    assert(all([os.path.exists(f) for f in fnames]))
    ms, names = zip(*sorted([(json.load(open(fname))['m'], fname) for fname in fnames if '_expansions' not in fname]))
    return names

def get_fnames_edge_expansions():
    Wedge_files = glob.glob(f'{root}/Wedges_pass*.json')
    def get_pass(f):
        match = re.search(r'Wedges_pass([0-9]+).json', f).group(1)
        return int(match)
    Wedge_filess = sorted(Wedge_files, key=get_pass)
    fname = Wedge_filess[-1]
    import os
    fnames = [f'{root}/PCSAFT_VLE_m{1/f:0.12e}_expansions.json' for f in json.load(open(fname))["Wedges"]]
    assert(all([os.path.exists(f) for f in fnames]))
    return fnames, json.load(open(fname))["Wedges"]

def plot_intervals():
    fig, ax = plt.subplots(1,1,figsize=(3.5, 2.5))
    for f, w in zip(*get_fnames_edge_expansions()):
        j = json.load(open(f))
        edges = [ex['xmin'] for ex in j]
        edges += [ex['xmax'] for ex in j]
        edges = list(set(sorted(edges)))
        plt.plot(w*np.ones_like(edges), edges, 'o')
    plt.gca().set(xlabel=r'$w=1/m$', ylabel=r'$\Theta$')
    plt.tight_layout(pad=0.2)
    plt.savefig('intervals.pdf')
    plt.show()

def plot_all_VLE():

    fig, ax = plt.subplots(1,1,figsize=(3.3,3))
    for f in get_fnames_edges():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        if m > 64 and m != 80:
            continue
        Ttildemin = np.min(df['Ttilde'])
        
        Tred = df['Ttilde']/j['Ttildec']
        line, = plt.plot(df['rhotildeL']/j['rhotildec'], Tred, label='')
        plt.plot(df['rhotildeV']/j['rhotildec'], Tred, color=line.get_color())
    plt.axvline(1.0, dashes=[1,1])
    plt.gca().set(
        xlabel=r'$\widetilde{\rho}/\widetilde{\rho}_{\rm crit}$', 
        ylabel=r'$\widetilde{T}/\widetilde{T}_{\rm crit}$'
    )
    # plt.legend(loc='best', ncol=3, fontsize=4)
    plt.tight_layout(pad=0.2)
    plt.savefig('PCSAFT_all_VLE.pdf')
    plt.close()

def plot_Tmin():
    ooo = []
    fig, ax = plt.subplots(1,1,figsize=(3.3,3))
    for f in get_fnames_edges():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        Tred = df['Ttilde']/j['Ttildec']
        line, = plt.plot(df['rhotildeL']/j['rhotildeV'], Tred, label='')
        try:
            Tmin_red = float(scipy.interpolate.interp1d(df['rhotildeL']/df['rhotildeV'], Tred)(1e20))
            ooo.append({'m': j['m'], 'Tmin_red': Tmin_red})
        except:
            pass
    # plt.legend(loc='best')
    plt.xscale('log')
    plt.xlim(1, 1e40)

    plt.gca().set(
        xlabel=r"$\widetilde{\rho}'/\widetilde{\rho}''$",
        ylabel=r'$\widetilde{T}/\widetilde{T}_{\rm crit}$'
    )
    plt.xticks(10.0**np.array([0,6,12,18,24,30,36]))
    plt.tight_layout(pad=0.2)
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
    fig, (axV, axL) = plt.subplots(1, 2, figsize=(3.3,3), sharey=True)
    for f in get_fnames_edges():
        j = json.load(open(f))
        arrays = ['Ttilde', 'rhotildeL', 'rhotildeV']
        df = pandas.DataFrame({k:j[k] for k in arrays})
        m = j['m']
        if m > 80:
            continue
        
        # From the empirical fit
        c = [ 0.37627892, -2.20078778 ]
        Tredmin = np.exp(c[1])*np.power(m, c[0])
        Ttildemin = Tredmin*j['Ttildec']

        Theta = (df['Ttilde']-Ttildemin)/(j['Ttildec']-Ttildemin)
        
        line, = axL.plot(df['rhotildeL']/j['rhotildec'], Theta, label='')
        axV.plot(df['rhotildeV']/j['rhotildec'], Theta, color=line.get_color())

    axV.set(
        ylabel=r'$\Theta=(\widetilde{T}-\widetilde{T}_{\rm min})/(\widetilde{T}_{\rm crit}-\widetilde{T}_{\rm min})$'
    )
    axV.set(
        xlabel=r"$\widetilde{\rho}''/\widetilde{\rho}''_{\rm crit}$"
    )
    axL.set(
        xlabel=r"$\widetilde{\rho}'/\widetilde{\rho}'_{\rm crit}$"
    )
    axV.set_xscale('log')
    axV.set_xlim(1e-20,1)
    axV.set_xticks(10.0**np.array([0, -10.0, -20]))
    plt.ylim(0, 1)
    axL.set_xlim(1, 7)
    plt.tight_layout(pad=0.2)
    plt.savefig('PCSAFT_normalized_VLE.pdf')
    plt.close()

def plot_allrhoerr(Nm, fitted=False):
    for domain_index in range(8):
        if fitted:
            df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_fitted_Nm16.json')))
        else:
            df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_worstcase_Nm{Nm}.json')))
        err = df['rhotildeL_VLEmp']/df['rhotildeL_anc']-1
        err[err==0] = 1e-16
        baddies = ~np.isfinite(err)
        if sum(baddies) > 0:
            print(sum(baddies))
        sc = plt.scatter(df['Theta'], np.log10(np.abs(err)), c=df['1/m'], vmin=1/64, vmax=1/1, cmap='viridis', lw=0)
    
    plt.gca().set(xlabel=r'$\Theta$', ylabel=r'$\log_{10}(|err|)$', xlim=(0,1))
    cb = plt.colorbar()
    cb.set_label(r'$1/m$')
    plt.tight_layout(pad=0.2)
    if fitted:
        plt.savefig(f'all_fitted_devplot.pdf')
    else:
        plt.savefig(f'all_Nm{Nm}_devplot.pdf')
    plt.show()

def plot_allrhoerr_m(root, Nm, fitted=False, Theta_cutoff=0.0, suffix=''):
    fig, axes = plt.subplots(1,2,figsize=(6, 5),sharey=True, sharex=True)
    maxerr = 0
    for domain_index in range(10):
        if fitted:
            df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_fitted_Nm16.json')))
        else:
            df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_worstcase_Nm{Nm}.json')))
        df = df[df['Theta'] >= Theta_cutoff]
        df = df[df['Theta'] <= 1.0-2.2e-14]

        keys = ['rhotildeL','rhotildeV']
        for key, ax in zip(keys, axes):
            err = df[f'{key}_VLEmp']/df[f'{key}_anc']-1
            err[err==0] = 1e-16
            maxerr = max(maxerr, np.max(err))
            baddies = ~np.isfinite(err)
            badTheta = df.loc[baddies,'Theta']
            if sum(baddies) > 0:
                print(sum(baddies), np.min(badTheta))
            sc = ax.scatter(df['1/m'], np.log10(np.abs(err)), c=df['Theta'], vmin=Theta_cutoff, vmax=1, cmap='viridis', lw=0)
            ax.set_title(key)
    
    axes[1].set(xlabel=r'$1/m$', xlim=(1/64,1))
    axes[0].set(xlabel=r'$1/m$', ylabel=r'$\log_{10}(|\rho^\alpha_{\rm mp}/\rho^\alpha_{\rm SA}-1|)$', xlim=(1/64,1))
    cb = plt.colorbar(sc, ax=axes[1])
    cb.set_label(r'$\Theta$')
    axes[1].set_ylim(bottom=np.log10(9e-17), top=np.log10(maxerr*1.01))
    plt.tight_layout(pad=0.2)
    if fitted:
        plt.savefig(f'all_fitted_devplot_funcm{suffix}.pdf')
    else:
        plt.savefig(f'all_Nm{Nm}_devplot_funcm{suffix}.pdf')
    plt.close()



def plot_rhoerr(root, domain_index, Nm):
    if Nm == 0:
        df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_fitted_Nm16.json')))
    else:
        df = pandas.DataFrame(json.load(open(root+f'/{domain_index}PCSAFT_VLE_check_worstcase_Nm{Nm}.json')))
    
    # sc = plt.scatter(df['Theta'], df['1/m'], c=np.log10(np.abs(df['rhotildeL_VLE']/df['rhotildeL_anc']-1)), vmin=-16, vmax=2)
    # plt.gca().set(xlabel=r'$\Theta$', ylabel=r'$1/m$')
    # cb = plt.colorbar()
    # cb.set_label(r'$\log10(|err|)$')
    # if Nm == 0:
    #     plt.savefig(f'fitted_Thetam.pdf')
    # else:
    #     plt.savefig(f'Nm{Nm}_Thetam.pdf')
    # plt.show()

    sc = plt.scatter(df['Theta'], np.log10(np.abs(df['rhotildeL_VLE']/df['rhotildeL_anc']-1)), c=df['1/m'])
    plt.gca().set(xlabel=r'$\Theta$', ylabel=r'$\log10(|relerr|)$')
    cb = plt.colorbar()
    cb.set_label(r'$1/m$')
    
    if Nm == 0:
        plt.savefig(f'fitted_devplot.pdf')
    else:
        plt.savefig(f'Nm{Nm}_devplot.pdf')
    plt.show()

if __name__ == '__main__':
    import ChebTools
    
    # ce = ChebTools.ChebyshevExpansion(np.linspace(0,1,17), 1/64, 1/1)
    # print(1/ce.get_nodes_realworld())
    # quit()

    # plot_critical_curve()
    # plot_critical_curvedev()
    # plot_all_VLE()
    # plot_normalized_VLE()
    # plot_intervals()
    # plot_Tmin()
    # plot_allrhoerr(16, fitted=True)
    # plot_allrhoerr(16, fitted=False)
    plot_allrhoerr_m('bld', 16, fitted=True)
    plot_allrhoerr_m('bld', 16, fitted=False)
    # plot_allrhoerr_m('bld', 16, fitted=True, Theta_cutoff=0.9999, suffix='_nearcrit')