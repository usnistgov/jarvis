from pymatgen.util.plotting_utils import get_publication_quality_plot
import matplotlib,yaml,os
from pymatgen.io.vasp.inputs import Potcar,Incar, Kpoints
from pymatgen.electronic_structure.plotter import BSPlotterProjected, BSPlotter, DosPlotter
import matplotlib
matplotlib.use('Agg')
from pymatgen.io.vasp.outputs import Vasprun


run = Vasprun(ru,occu_tol=1e-2)
complete_dos = run.complete_dos
totup=complete_dos.get_densities(spin=Spin.up)
totdn=complete_dos.get_densities(spin=Spin.down)
en = complete_dos.energies
bandgap=str(run.eigenvalue_band_properties[0])
ef = complete_dos.efermi
en[:] = [x - ef for x in en]



#allpts = []
#x=complete_dos.densities
#y=complete_dos.energies

plt = get_publication_quality_plot(14, 10)
plt.xlabel("Energies (eV)")
plt.ylabel("DOS (arb. unit.)")
plt.plot(en,totup,'b',linewidth=3)
plt.plot(en,-totdn,'r',linewidth=3)
filename=str(storedir)+str("/")+str("TDos.png")
plt.xlim(-5,10)
plt.tight_layout()
plt.savefig(filename)
plt.close()

