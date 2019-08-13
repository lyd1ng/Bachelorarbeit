Dieses Verzeichnis beinhaltet
einen Test mit dem geprueft werden
sollte, ob float-Vektoren uber ein
Buffer an einen Kernel uebergeben werden
koennen. Das ist der Fall, jedoch muss
hostseitig immer ein float4-Vektor
verwendet werden, auch wenn auf der Seite
des Kernels ein Vektor geringerer Ordnung
definiert ist.
Weiterhin ist darauf zu achten, dass kein
bytepadding stattfindet (#pragma packed)
