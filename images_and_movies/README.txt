Dieses Verzeichnis zeigt ein paar Plots, die im Entwicklungsprozess der
FDTD-Implementierung in Python entstanden:

1d_pml_efactor_plot.png -- Der pml-Vorfaktor fuer das E-Feld
1d_pml_hfactor_plot.png -- Der pml-Vorfaktor fuer das H-Feld
constant_sin_instability.mp4 -- Eine ungleiche Ausbreitung ist fuer einen konstanten Sinus zu erkennen, dessen Spektrum illegale Frequenzen beinhaltet
waveguide_wrong.mp4 -- Ein erster Versuch einen Hohlleiter zu simulieren, der jedoch nicht korrekt ist. Der linke und obere Rand des Hohlleiters sind durch
		       einen perfekten magnetischen Leiter anstatt eines perfekten magnetischen Leiters begrenzt. Das folgt direkt aus der 
		       Dirichlet-Randbedinung.
waveguide_correct_efield.mp4 -- Eine korrekte Simulation eines Hohlleiters indem ein perfekter elektrischer Leiter an den oberen und linken Rand gesetzt wurde (elektrisches Feld)
waveguide_correct_pfield.mp4 -- Eine korrekte Simulation eines Hohlleiters indem ein perfekter elektrischer Leiter an den oberen und linken Rand gesetzt wurde (poynting Feld)

