/*****************************************************************************
 *			         XY Model 2D			             *
 *			       Cone Interactions			     *
 *			        Pedro H Mendes 			             *
 ****************************************************************************/
/*****************************************************************************
 *	Compile:							     *
 *	gcc -Wall xy_cone.c -lm -O3	 			             *
 *									     *
 *	Flags:								     *
 *	-DDATA -> temporal series of energy				     *
 *	-DGNU  -> gnuplot visualization					     *
 *									     *
 *	Run:								     *
 *	./a.out TEMP DEG SEED						     *
 *	./a.out TEMP DEG SEED | gnuplot 				     *
 *****************************************************************************/

/*****************************************************************************
 *                             	   INCLUDES                                  *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"mc.h"

/*****************************************************************************
 *                               DEFINITIONS                                 *
 ****************************************************************************/
#define 		L			32
#define 		L2 	 		(L*L)
#define 		TRAN			0
#define 		TMAX			100000
#define			N			4.0
#define 		MIN(a,b) 		(((a)<(b))?(a):(b))

/*****************************************************************************
 *                           GLOBAL VARIABLES                                *
 ****************************************************************************/
double dE, M, ET, U, CV;
double EE, ET2, M2, M4;
int J, THETA_DEGREES;
double THETA;

/*****************************************************************************
 *                              FUNCTIONS                                    *
 ****************************************************************************/
void initialize(double *spin, int **neigh);
void mc_routine(double *spin, int **neigh, double TEMP);
void print_states(double *spin, double TEMP, int choice);
void gnuplot_view(int tempo, double *spin);
void calculate_quantities(double *spin, double TEMP);

/*****************************************************************************
 *                             MAIN PROGRAM                                  *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	clock_t t_i, t_f;
	t_i = clock();

	int mcs;
	int i;
	double TEMP, CPU_TIME;

	TEMP = atof(argv[1]) / 100;
	THETA_DEGREES = atof(argv[2]);
	THETA = 2*M_PI*(1.0*THETA_DEGREES/360);

	double *spin;
	int **neigh;

	spin = (double*)malloc(L2*sizeof(double));
	neigh = (int**)malloc(L2*sizeof(int*));

	for(i=0; i<L2; i++)
        {
                neigh[i] = (int*)malloc(4*sizeof(int));
        }

	seed = start_randomic_seed_ext(atoi(argv[3]));
	initialize(spin, neigh);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		mc_routine(spin, neigh, TEMP);

	}

#ifdef DATA
	char Arq1[100];
	FILE *arq1;

	sprintf(Arq1, "temp_T%.3lfTheta=%.0lfL%dS%ld.dat", TEMP, atof(argv[2]), L, seed);
	arq1 = fopen(Arq1, "w");
	//fprintf(arq1, "#seed = %ld\n#MCS,ET,M,M2,M4,U,ET2,CV,EE\n", seed);
	fprintf(arq1, "#seed = %ld\n#MCS,M,U\n", seed);
#endif

	for(mcs=0; mcs<TMAX; mcs++)
	{
		mc_routine(spin, neigh, TEMP);
#ifdef GNU
        gnuplot_view(mcs,spin);
#endif

#ifdef DATA
	//fprintf(arq1, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n", mcs, ET, M, M2, M4, U, ET2, CV, EE);
	fprintf(arq1, "%d,%f,%f\n", mcs, M, U);
#endif
	}
	EE += ET*ET;

	print_states(spin, TEMP, 0);

#ifdef DATA
	fclose(arq1);
#endif

	free(spin);

	for(i=0; i<L2; i++)
	{
		free(neigh[i]);
	}

	t_f = clock();
	CPU_TIME = (double)(t_f - t_i)/CLOCKS_PER_SEC;

	printf("#%lf\n", CPU_TIME);

	return 0;
}

/*****************************************************************************
 *                             INITIALIZATION                                *
 ****************************************************************************/
void initialize(double *spin, int **neigh)
{
	int i;

	for(i=0; i<L2; i++)
	{
		spin[i] = 2*M_PI*FRANDOM;
	}

	for(i=0; i<L2; i++)
	{
		neigh[i][1] = (i-L+L2)%L2;              //up
        neigh[i][0] = (i+1)%L + (i/L)*L;        //right
        neigh[i][3] = (i+L)%L2;                 //down
        neigh[i][2] = (i-1+L)%L + (i/L)*L;      //left
	}

	return;
}

/*****************************************************************************
 *                     	      CALCULATE QUANTITIES                           *
 ****************************************************************************/
void calculate_quantities(double *spin, double TEMP){
    // eu sou jesus
    double sum_sines = 0;
    double sum_cosines = 0;
    double sum_sines_squared = 0;
    double sum_cosines_squared = 0;
    double sum_sines_fourth = 0;
    double sum_cosines_fourth = 0;
    int    i;

    for(i = 0; i < L2; i++){
      sum_sines           += sin(spin[i]);
      sum_cosines         += cos(spin[i]);

      sum_sines_squared   += sin(spin[i]) * sin(spin[i]);
      sum_cosines_squared += cos(spin[i]) * cos(spin[i]);

      sum_sines_fourth    += sin(spin[i]) * sin(spin[i]) * sin(spin[i]) * sin(spin[i]);
      sum_cosines_fourth  += cos(spin[i]) * cos(spin[i]) * cos(spin[i]) * cos(spin[i]);
    }

    M  = sqrt(sum_sines*sum_sines + sum_cosines * sum_cosines) / L2;
    M2 = sqrt(sum_sines_squared*sum_sines_squared + sum_cosines_squared*sum_cosines_squared) / (L2);
    M4 = sqrt(sum_sines_fourth*sum_sines_fourth + sum_cosines_fourth*sum_cosines_fourth) / (L2);

    //binder cummulant
    U  = 1 - (M4)/(3 * M2*M2);

    CV = fabs((ET2 / L2) - (ET / L2)*(ET / L2)) / (TEMP*TEMP);

    //printf("mean energy squared per spin = %f, (mean energy) squared per spin = %f \n", (ET2 / L2), (ET / L2)*(ET / L2));
}

/*****************************************************************************
 *                     	      MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void mc_routine(double *spin, int **neigh, double TEMP)
{
	int i, j, t, total_J, int_flip;
	double Ei, Ef, G, final_individual_energy;
	double vi, D_ang, min_ang;
	double flip;

	for(t=0; t<L2; t++)
	{
		i = FRANDOM*L2;

		Ei = 0;
		total_J  = 0;
		for(j=0; j<4; j++)
		{
            vi = ((j)*2*M_PI)/4;
            //printf("spin_i = %lf ", spin[i]);
            //printf("-->vi_%d = %lf\n", j,vi);
            D_ang = fabs(spin[i]-vi);
            min_ang = MIN((2*M_PI) - D_ang,D_ang);
            //printf("min_ang --> %lf\n", min_ang);
            //otimizaçao da parte de baixo
            //J = ((min_ang) < THETA/2);
            //printf("theta/2 = %lf\n", THETA/2);

            if(min_ang < THETA/2)
            {
                    J=1;
            }
            else
            {
                    J=0;
            }
			//printf("J --> %d\n", J);
			total_J += J;
                        Ei += J*(cos(spin[i]-spin[neigh[i][j]]));
		}

		if (total_J > 0) {
		  //printf("%d \n", total_J);
		}

		flip = (spin[i] + (2*M_PI)*FRANDOM);

		int_flip = (360 * (flip / (2*M_PI)));

		int_flip = int_flip % 360;

		flip     = (int_flip * 1.0 * (2*M_PI/360));

		//printf("%f \n", flip);

		Ef = 0;
		for(j=0; j<4; j++)
		{
			vi = 1+j;
                        vi = (vi*2*M_PI)/4;
                        D_ang = fabs(flip-vi);
                        min_ang = MIN((2*M_PI) -D_ang,D_ang);

                        //otimizaçao da parte de baixo
                        J = (min_ang < THETA/2);

                            /*if(min_ang < THETA/2)
                        {
                                J=1;
                        }
                        else
                        {
                                J=0;
                        }*/

                        Ef += J*(cos(flip-spin[neigh[i][j]]));
		}

		dE = Ei-Ef;

		G = (1-tanh(dE/(2*TEMP)))/2;

		if(FRANDOM < G)
		{
			spin[i] = flip;
		}
	}

	ET  = 0;
	ET2 = 0;

	for(i=0; i<L2; i++)
	{
	    final_individual_energy = 0;
		for(j=0;j<4;j++)
		{
			vi = 1+j;
			vi = (vi*2*M_PI)/4;
			D_ang = fabs(spin[i]-vi);
			min_ang = MIN((2*M_PI) -D_ang,D_ang);

			//otimizaçao
            J = (min_ang < THETA/2);

			final_individual_energy += J * (cos(spin[i]-spin[neigh[i][j]]));
		}
		ET  += final_individual_energy;
		ET2 += final_individual_energy * final_individual_energy;
	}

	EE = ET * ET;

	calculate_quantities(spin, TEMP);

	return;
}

/*****************************************************************************
 *                                PRINT STATES                               *
 ****************************************************************************/
void print_states(double *spin, double TEMP, int choice)
{
        int i;

        if(choice == 0)
                {
                        char fp[100];
                        FILE *fp1;
		sprintf(fp, "finalconfig_T%.3lfL%dA%dS%ld.dat", TEMP, L, THETA_DEGREES, seed);
                        fp1 = fopen(fp, "w");

		for(i=0; i<L2; i++)
		{
			fprintf(fp1, "%f ", spin[i]);
		}

        /*
		esse jeito printa em formato LxL quadrado
                        for(i=0; i<L; i++)
                        {
			for(j=0; j<L; j++)
                                {
                                        fprintf(fp1, "%f ", spin[i + j*L]);
                                }
                                fprintf(fp1, "\n");
                        }
        */

                        fclose(fp1);
                }

               	return;
}

/**************************************************************************
 *                       Visualization  routine                           *
 *      compile with: -DGNUPLOT                                           *
 *      use as: ./a.out  *kwargs| gnuplot                                 *
 *************************************************************************/
void gnuplot_view(int tempo, double *spin)
{
        int i,j,sitio;

        printf("set title \'tempo: %d \' \n",tempo);
        printf("set size square\n");
        printf("plot \'-\' matrix with image\n");

        for(j=0;j<L;j++)
        {
                for(i=0;i<L;i++)
                {
                        sitio = i+j*L;
                        printf(" %f",spin[sitio]);
                }
                printf("\n");
        }

        printf("e\n pause 0.05\n");
        printf("\n\n");

        return;
}
