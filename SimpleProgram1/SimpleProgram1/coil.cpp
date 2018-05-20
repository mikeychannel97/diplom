#include "stdafx.h"


/* массив для отладочных выводов в частных программах  */
int    M_PECH[50];
/* массив для отладочных выводов в inskor  */
int mpint[40];
/* массив для отладочных выводов в INTU_KOR  */
int M_P_INTK[10];


double MKC[35];/* коор. 7 постов */
double TKOC[200];/* константы */

void coil(void)
	{
	 int Kpost;
	 int N_TKOC, N_TKOCmax;
	 int i,j;
	 char  kolich [100];
     FILE *in;

/*  чтение MKC из файлa  MKC.txt*/
	 fopen_s(&in, "MKC.txt", "r");
     /*
	 if ((in = fopen("MKC.txt", "r")) == NULL)
	 {
	     fprintf(stderr, "1 Cannot open input file MKC.txt.txt.\n");
	     exit(0);		
	 }
     */
	 fgets (kolich,100,in);
	 fscanf_s(in, "%d \n", &Kpost);                            /* Kpost */
     fgets (kolich,100,in);                           
     for (i=0; i<Kpost*5;i++)
	 {
         fscanf_s(in, " %s \n", kolich, 100); 
/* 03.10.2006 */
         MKC[i]   = atof ( kolich);
	 }
     fclose(in);/* MKC.txt  */

/*  чтение TKOC[] из файлa  TKOC_2000.txt*/
	 fopen_s(&in, "TKOC_2000.txt", "r");
     /*
	 if ((in = fopen("TKOC_2000.txt", "r")) == NULL)
	 {
	     fprintf(stderr, "1 Cannot open input file TKOC_2000.txt.txt.\n");
	     exit(0);		
	 }
     */

     fgets (kolich,100,in);
     fscanf_s(in, "%d \n", &N_TKOCmax);                        /* N_TKOCmax */
		for (i=0; i<N_TKOCmax; i++)
		{
           TKOC[i] = 1.0e+17;
		}
     fgets (kolich,100,in);
     fscanf_s(in, "%d \n", &N_TKOC);                           /* N_TKOC */
     for (i=0; i<N_TKOC;i++)
	 {
         fgets (kolich,100,in);
		 fscanf_s(in, " %d  %s \n", &j, kolich, 100); 
         TKOC[i]   = atof ( kolich);
	 }
     fclose(in);/* TKOC_2000.txt  */


/* обнуление массива признаков печати для программы int92 */
	 	for (i=0; i<40;i++)
		{
		     mpint[i] = 0;
		}
/* обнуление массивa признаков печати для программы INTU_KOR  */
	 	for (i=0; i<10;i++)
		{
		     M_P_INTK[i] = 0;
		}
	 	for (i=0; i<50;i++)
		{
		     M_PECH[i] = 0;
		}
	   return;
}