Allgemeinse:
============
Dieses Verzeichnis beinhaltet alles was im Rahmen meiner
Bachelorarbeit zustande kam. Das sind deutlich mehr als
tatsaechlich in die Bachelorthesis kamen, da ich die
FDTD fuer den 3d-Fall schrittweise implementieren wollte
sich jedoch nur die Implementierung fuer den 3d-Fall nutzen
laesst um Drahtantennen zu simulieren.

Inhalt:
=======
benchmarking_results: Beinhaltet die Benchmarking Ergebnisse python zu c und c zu opencl
fdtd_python: Erste Prototypen der FDTD in Python
opencl: Beinhaltet die Implementierungen der FDTD in OpenCL sowie alle anderen Programme die ich schrieb
	um mit OpenCL vertraut zu werden
reference_implementation: Die C Referenzimplementationen der FDTD
thesis: Beinhaltet Bachelorarbeit als pdf.

Bilder und Videos:
==================
Da die Python-Implementierungen schoene Bilder und Videos machten
existiert ein Ordner "images_and_movies".

Ausfuehren der Python Programme:
================================
Es wird das Modul matplotlib benoetigt.
Ansonsten sollten keine weiteren Abhaengigkeiten vorhanden sein.

Compilieren der C Programme:
============================
Zu jedem C-Programm steht eine Makefile bereit,
es wird jedoch mitunter die OpenCL Bibliothek benoetigt.
