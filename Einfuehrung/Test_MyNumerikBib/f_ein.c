#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <MyNumerikBib.h>
/* Einlesen von Daten */
void read_array(int n,double *b)
{
	FILE   *fd;		/* Filedeskriptor */ 
	int i;
        /* Einlesen von Daten */
	fd = fopen ("Vector.dat", "r");
	//Pruefen, ob die Datei geoeffnet werden kann
	if (fd == NULL)
	{
	    fprintf (stderr, "Datei Vector.dat konnte nicht zum Lesen geöffnet werden!\n");
	    exit(1);
	}
	/* Einlesen des Vektors*/
	for (i=0; i<n; i++){
		fscanf(fd, "%lf*[^\n]", &b[i]);
	}
	/* Filedeskriptor schliessen */
	fclose(fd);
	
}
