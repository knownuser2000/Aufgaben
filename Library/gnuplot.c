// Created by Stefan Jakowinow

#include "gnuplot.h"

void gnuplot(struct gnuplot_arg arg, char **legend, double *x, int num, ...)
{
    // ______Failsafe_______
    // checking if enough arrays are given
    if(num < 1){printf("\t ERROR: Not enough arrays!\n"); exit(1);}
    // check if length is given
    if(arg.length <= 0){printf("\t ERROR: False Length!\n"); exit(1);}
    // create standard filename if not provided
    if (arg.inputFile == NULL){arg.inputFile = "inputGnuplot.dat";}

    double *y,*arr;
    arr = calloc(arg.length,sizeof(double));
    y = calloc(arg.length*num,sizeof(double));

    va_list list;
    va_start(list, num);

    for(int i=0;i<num;i++){
        arr = va_arg(list,double *);
        for (int j = 0; j < arg.length; j++){
            y[j*num +i] = arr[j];
        }
    }

    va_end(list);

    createInputFile(arg, x, num, y);

    createGnuplotFile(arg, legend, num);

    printf("\t Creating Gnuplot!\n");
    system("gnuplot gnuplot.gpl -p");
    //printf("\t Gnuplot created!\n\n");

    //free(arr); free(x); free(y);
}

void createInputFile(struct gnuplot_arg arg, double *x, int num, double *y)
{
    printf("\n\t Creating input file for gnuplot!\n");

    FILE   *fd;
	fd = fopen (arg.inputFile, "w");

	if (fd == NULL){
	    fprintf (stderr, "Datei %s konnte nicht zum Schreiben geöffnet werden!\n", arg.inputFile);
	    exit(1);
	}

	for(int i = 0; i<arg.length ; i++){
        fprintf(fd,"%E\t",x[i]);
        for(int j = 0; j<num ; j++){
            fprintf(fd,"%E\t",y[i*num+j]);
        }
        fprintf(fd,"\n");
	}

	fclose(fd);

    //printf("\t Input file created!\n\n");
}

void createGnuplotFile(struct gnuplot_arg arg, char **legend, int num)
{
    printf("\t Creating gnuplot file!\n");

    FILE *fd;
    char fileName[] = "gnuplot.gpl";
    fd = fopen(fileName, "w");
    if (fd == NULL){
        fprintf(stderr, "Datei %s konnte nicht zum Schreiben (GNu) geöffnet werden!\n", fileName);
        exit(1);
    }

    if (arg.grid == 1){
        fprintf(fd, "set grid\n");}
    fprintf(fd, "set zeroaxis\n");
    fprintf(fd, "set title \"%s\"\n", arg.plotTitle != NULL ? arg.plotTitle : "Gnuplot Title");
    fprintf(fd, "set xlabel \"%s\"\n", arg.xlabel != NULL ? arg.xlabel : "x Axis");
    fprintf(fd, "set ylabel \"%s\"\n", arg.ylabel != NULL ? arg.ylabel : "y Axis");
    fprintf(fd, "set style data %s \n", arg.linestyle != NULL ? arg.linestyle : "lines");
    fprintf(fd, "set datafile separator \"\\t\" \n");
    fprintf(fd, "set zeroaxis\n");

    fprintf(fd, "plot \"%s\" using 1:2 title \"%s\"\\\n", arg.inputFile, 
            legend[0] != NULL ? legend[0] : " ");

    for(int i=0;i<num-1;i++){
        fprintf(fd, "\t\t, '' using 1:%i title \"%s\"\\\n", i+3, legend[i+1] != NULL ? legend[i+1] : " ");
    }

    fclose(fd);
    //printf("\t Gnuplot file created!\n\n");
}

// matrix variant
void gnuplot_matrix(struct gnuplot_arg arg, char **legend, double *x, int num, double *y_matrix)
{
    // ______Failsafe_______
    // checking if enough arrays are given
    if(num < 1){printf("\t ERROR: Not enough arrays!\n"); exit(1);}
    // check if length is given
    if(arg.length <= 0){printf("\t ERROR: False Length!\n"); exit(1);}
    // create standard filename if not provided
    if (arg.inputFile == NULL){arg.inputFile = "inputGnuplot.dat";}


    createInputFile(arg, x, num, y_matrix);

    createGnuplotFile(arg, legend, num);

    printf("\t Creating Gnuplot!\n");
    system("gnuplot gnuplot.gpl -p");
}

// plot with 2 different y axis
void gnuplot_multi_y(struct gnuplot_arg arg, char **legend_y1, char **legend_y2,
                     double *x, int num_y1, int num_y2, ...)
{
    // ______Failsafe_______
    // checking if enough arrays are given
    if(num_y1 <= 0 || num_y2 <= 0){printf("\t ERROR: Not enough arrays!\n"); exit(1);}
    // check if length is given
    if(arg.length <= 0){printf("\t ERROR: False Length!\n"); exit(1);}
    // create standard filename if not provided
    if (arg.inputFile == NULL){arg.inputFile = "inputGnuplot.dat";}

    double *y1, *y2, *arr;
    arr = calloc(arg.length,sizeof(double));
    y1 = calloc(arg.length*num_y1,sizeof(double));
    y2 = calloc(arg.length*num_y2,sizeof(double));

    va_list list;
    va_start(list, num_y2);

    // saving for y1
    for(int i=0;i<num_y1;i++){
        arr = va_arg(list,double *);
        for (int j = 0; j < arg.length; j++){
            y1[j*num_y1 +i] = arr[j];
        }
    }
    // saving for y2
    for(int i=0;i<num_y2;i++){
        arr = va_arg(list,double *);
        for (int j = 0; j < arg.length; j++){
            y2[j*num_y2 +i] = arr[j];
        }
    }
    
    va_end(list);

    createInputFile_multi_y(arg, num_y1, num_y2, x, y1, y2);

    createGnuplotFile_multi_y(arg, legend_y1, legend_y2, num_y1, num_y2);

    printf("\t Creating Gnuplot!\n");
    system("gnuplot gnuplot_multi_y.gpl -p");
    //printf("\t Gnuplot created!\n\n");

    //free(arr); free(x); free(y);
}

void createInputFile_multi_y(struct gnuplot_arg arg, int num_y1, int num_y2, double *x, double *y1, double *y2)
{
    printf("\n\t Creating input file for gnuplot!\n");

    FILE   *fd;
	fd = fopen (arg.inputFile, "w");

	if (fd == NULL){
	    fprintf (stderr, "Datei %s konnte nicht zum Schreiben geöffnet werden!\n", arg.inputFile);
	    exit(1);
	}

    for (int i = 0; i < arg.length; i++){
        // x axis
        fprintf(fd, "%E\t", x[i]);
        // y1 axis
        for (int j = 0; j < num_y1; j++){
            fprintf(fd, "%E\t", y1[i * num_y1 + j]);
        }
        // y2 axis
        for (int j = 0; j < num_y2; j++){
            fprintf(fd, "%E\t", y2[i * num_y2 + j]);
        }
        fprintf(fd, "\n");
    }
	fclose(fd);
    //printf("\t Input file created!\n\n");
}

void createGnuplotFile_multi_y(struct gnuplot_arg arg, char **legend_y1, char **legend_y2, int num_y1, int num_y2)
{
    printf("\t Creating gnuplot file!\n");

    FILE *fd;
    char fileName[] = "gnuplot_multi_y.gpl";
    fd = fopen(fileName, "w");
    if (fd == NULL){
        fprintf(stderr, "Datei %s konnte nicht zum Schreiben (GNu) geöffnet werden!\n", fileName);
        exit(1);
    }

    if (arg.grid == 1){
        fprintf(fd, "set grid\n");}
    fprintf(fd, "set zeroaxis\n");
    fprintf(fd, "set title \"%s\"\n", arg.plotTitle != NULL ? arg.plotTitle : "Gnuplot Title");
    fprintf(fd, "set xlabel \"%s\"\n", arg.xlabel != NULL ? arg.xlabel : "x Axis");
    fprintf(fd, "set ylabel \"%s\"\n", arg.ylabel != NULL ? arg.ylabel : "y1 Axis");
    fprintf(fd, "set y2tics \n");
    fprintf(fd, "set ytics nomirror \n");
    fprintf(fd, "set y2label \"%s\"\n", arg.y2label != NULL ? arg.y2label : "y2 Axis");
    fprintf(fd, "set style data %s \n", arg.linestyle != NULL ? arg.linestyle : "lines");
    fprintf(fd, "set datafile separator \"\\t\" \n");

    // 1 plot
    fprintf(fd, "plot \"%s\" using 1:2 axis x1y1 title \"%s\"\\\n", arg.inputFile, legend_y1[0] != NULL ? legend_y1[0] : " ");

    // every other y1 plot
    for(int i=0; i < num_y1-1; i++)
    {
        fprintf(fd, "\t\t, '' using 1:%i axis x1y1 title \"%s\"\\\n", i+3, legend_y1[i+1] != NULL ? legend_y1[i+1] : " ");
    }
    // y2 plot
    for(int i=0; i < num_y2; i++)
    {
        fprintf(fd, "\t\t, '' using 1:%i axis x1y2 title \"%s\"\\\n", i+num_y1+2, legend_y2[i] != NULL ? legend_y2[i] : " ");
    }

    fclose(fd);
    //printf("\t Gnuplot file created!\n\n");
}

void gnuplot_multi_2D_in_3d(struct gnuplot_arg arg, double *x, double *y, int num, double *z_matrix)
{
    // ______Failsafe_______
    // checking if enough arrays are given
    if(num < 1){printf("\t ERROR: Not enough arrays!\n"); exit(1);}
    // check if length is given
    if(arg.length <= 0){printf("\t ERROR: False Length!\n"); exit(1);}
    // create standard filename if not provided
    if (arg.inputFile == NULL){arg.inputFile = "inputGnuplot.dat";}

    createInputFile_multi_2D_in_3d(arg, x, num, z_matrix);

    createGnuplotFile_multi_2D_in_3d(arg, y, num);

    printf("\t Creating Gnuplot!\n");
    system("gnuplot \"gnuplot_multi_2D_in_3D.gpl\" -persist");
}

void createInputFile_multi_2D_in_3d(struct gnuplot_arg arg, double *x, int num, double *z_matrix)
{
    printf("\n\t Creating input file for gnuplot!\n");

    FILE   *fd;
	fd = fopen (arg.inputFile, "w");

	if (fd == NULL){
	    fprintf (stderr, "Datei %s konnte nicht zum Schreiben geöffnet werden!\n", arg.inputFile);
	    exit(1);
	}

	for(int i = 0; i<arg.length ; i++){
        fprintf(fd,"%E\t",x[i]);
        for(int j = 0; j<num ; j++){
            fprintf(fd,"%E\t",z_matrix[i*num+j]);
        }
        fprintf(fd,"\n");
	}

	fclose(fd);
}

void createGnuplotFile_multi_2D_in_3d(struct gnuplot_arg arg, double *y, int num)
{
    printf("\t Creating gnuplot file!\n");

    FILE *fd;
    char fileName[] = "gnuplot_multi_2D_in_3D.gpl";
    fd = fopen(fileName, "w");
    if (fd == NULL){
        fprintf(stderr, "Datei %s konnte nicht zum Schreiben (GNu) geöffnet werden!\n", fileName);
        exit(1);
    }

    if (arg.grid == 1){
        fprintf(fd, "set grid\n");}
    fprintf(fd, "set zeroaxis\n");
    fprintf(fd, "set title \"%s\"\n", arg.plotTitle != NULL ? arg.plotTitle : "Gnuplot Title");
    fprintf(fd, "set xlabel \"%s\"\n", arg.xlabel != NULL ? arg.xlabel : "x Axis");
    fprintf(fd, "set ylabel \"%s\"\n", arg.ylabel != NULL ? arg.ylabel : "y Axis");
    fprintf(fd, "set zlabel \"%s\"\n", arg.zlabel != NULL ? arg.zlabel : "z Axis");
    fprintf(fd, "set key noautotitle\n");
    fprintf(fd, "set style data %s \n", arg.linestyle != NULL ? arg.linestyle : "lines");
    fprintf(fd, "set datafile separator \"\\t\" \n");
    fprintf(fd, "set zeroaxis\n");
    fprintf(fd, "set view %i,%i\n", arg.view_rot_x != 0 ? arg.view_rot_x : 50, arg.view_rot_z != 0 ? arg.view_rot_z : 30);

    fprintf(fd, "splot \"%s\" using 1:(%lf):2 \\\n", arg.inputFile, y[0]);

    for(int i=0;i<num-1;i++){
        fprintf(fd, "\t\t, '' using 1:(%lf):%i \\\n", y[i+1] ,i+3);
    }

    fclose(fd);
}