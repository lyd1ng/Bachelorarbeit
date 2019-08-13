# fdtd
A bunch of fdtd related programmes

fdtd_1d: Maxwellgleichungen 1d geloest. Mit perfekter absorbierender Randbedingung (das ist nur 1d moeglich)
fdtd_1d_pml: Maxwellgleichungen 1d geleost. Mit PML.

fdtd_2d: Maxwellgleichungen 2d geloest. Ohne absorbierende Randbedingung.
fdtd_2d_pml: Maxwellgleichungen 2d geloest. Mit PML

fft_ramped_sin: Eine numerische Spektralanayse eines ramped sinus, da ich sie analytisch nicht herleiten konnte.
		Leider ist die numerische Spektralanalyse nichtsaussagend

images_movies: Ein paar plots.


# PML (perfectly matched layer):
Das PML ist ein theoretisches anisotropes Material, dass Wellen aller Frequenzen und unter allen Winkeln reflektionsfrei
die Mediumsgrenze passieren laesst und Energie aus der Welle absorbiert.
Damit kann alle Energie aus einer Welle absorbiet werden bevor sie am Rande des Universums reflektiert wird.
Somit werden abgestrahlte wellen nicht in den sog. Problemraum zurueckreflektiert und der Freiraum simuliert.

Universum:

+-------------------------------------------------+
|                                               <-------------------PML
|   +-----------------------------------------+   |
|   |                                         |   |
|   |                                       <-----------------------Problem Space
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   |                                         |   |
|   +-----------------------------------------+   |
|                                                 |
+-------------------------------------------------+
