#pragma once

// Created by Stefan Jakowinow

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>


///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////2D-Plot////////////////////////////////////

// always initialize the struct as = {.grid = 1, ....}
struct gnuplot_arg{
    int grid;           // 0:: no grid; 1:: set grid
    int length;         // length of given arrays
    int view_rot_x;     // rotation angles (in degrees): x-axis
    int view_rot_z;     // rotation angles (in degrees): z-axis
    char *plotTitle;    // title of the plot
    char *xlabel;       // label for x axis
    char *ylabel;       // label for y axis
    char *y2label;      // label for 2. y axis
    char *zlabel;       // label for the z axis
    char *inputFile;    // name of a file
    char *linestyle;    // style of plotpoints
};

/*  2D-Plot can take any amount of arrays as input
    example:    gnuplot(plotstruct, legend, x, 3, cA, cB, cC)
    with legend:: char *legend[] = {"cA","cB","cC"}
         x:: array for x axis
         num:: number of arrays for y axis
         (...):: arrays for y axis
*/
void gnuplot(struct gnuplot_arg arg, char **legend, double *x, int num, ...);

void createInputFile(struct gnuplot_arg arg, double *x, int num, double *y);

void createGnuplotFile(struct gnuplot_arg arg, char **legend, int num);

void gnuplot_matrix(struct gnuplot_arg arg, char **legend, double *x, int num, double *y_matrix);

/* plot with 2 different y axis and same x axis
    example: gnuplot_multi_y(plotstruct, legend_y1, legend_y2, x, 2, 2, eps1, eps2, abs1, abs2);
    with x:: array for x axis
         num_y1 :: number of columns for y1 axis
         num_y2 :: number of columns for y2 axis
         (...):: num_y1 times arrays for y1 axis then num_y2 times arrays for y2 axis
*/
void gnuplot_multi_y(struct gnuplot_arg arg, char **legend_y1, char **legend_y2, double *x, int num_y1, int num_y2, ...);

void createInputFile_multi_y(struct gnuplot_arg arg, int num_y1, int num_y2, double *x, double *y1, double *y2);

void createGnuplotFile_multi_y(struct gnuplot_arg arg, char **legend_y1, char **legend_y2, int num_y1, int num_y2);

// plotting in 3D with x,y,z axis
void gnuplot_multi_2D_in_3d(struct gnuplot_arg arg, double *x, double *y, int num, double *z_matrix);

void createInputFile_multi_2D_in_3d(struct gnuplot_arg arg, double *x, int num, double *z_matrix);

void createGnuplotFile_multi_2D_in_3d(struct gnuplot_arg arg, double *y, int num);
