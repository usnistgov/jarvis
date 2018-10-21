from pymatgen.util.plotting_utils import get_publication_quality_plot
import matplotlib,yaml,os
from pymatgen.io.vasp.inputs import Potcar,Incar, Kpoints
from pymatgen.electronic_structure.plotter import BSPlotterProjected, BSPlotter, DosPlotter
import matplotlib,time
matplotlib.use('Agg')
from pymatgen.io.vasp.outputs import Vasprun



def banddos(pref='',storedir=None):
    ru=str("vasprun.xml")
    kpfile=str("KPOINTS")




    run = Vasprun(ru, parse_projected_eigen = True)
    bands = run.get_band_structure(kpfile, line_mode = True, efermi = run.efermi)
    bsp =  BSPlotter(bands)
    zero_to_efermi=True
    bandgap=str(round(bands.get_band_gap()['energy'],3))
    print "bg=",bandgap
    data=bsp.bs_plot_data(zero_to_efermi)
    plt = get_publication_quality_plot(12, 8)
    band_linewidth = 3
    x_max = data['distances'][-1][-1]
    print (x_max)
    for d in range(len(data['distances'])):
       for i in range(bsp._nb_bands):
          plt.plot(data['distances'][d],
                 [data['energy'][d]['1'][i][j]
                  for j in range(len(data['distances'][d]))], 'b-',
                 linewidth=band_linewidth)
          if bsp._bs.is_spin_polarized:
             plt.plot(data['distances'][d],
                     [data['energy'][d]['-1'][i][j]
                      for j in range(len(data['distances'][d]))],
                     'r--', linewidth=band_linewidth)
    bsp._maketicks(plt)
    if bsp._bs.is_metal():
         e_min = -10
         e_max = 10
         band_linewidth = 3

    for cbm in data['cbm']:
            plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                        s=100)

            for vbm in data['vbm']:
                plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                            s=100)


    plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
    ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi \
       else r'$\mathrm{Energy\ (eV)}$'
    plt.ylabel(ylabel, fontsize=30)
    plt.ylim(-4,4)
    plt.xlim(0,x_max)
    plt.tight_layout()
    plt.savefig('BAND.png',img_format="png")

    plt.close()






banddos()

