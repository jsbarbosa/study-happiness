#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 100000;
double dx = 0.0001, E, omega;
double *psi;
double psi_prime;

double *solver();
int main(int argc, char **argv)
{
    double acceptance = 0.5, current = 100;
    psi = malloc(N*sizeof(double));
    int i = 0, j = 0;
    
    FILE *functions = fopen("functions.dat","w");
    FILE *energies = fopen("energies.dat", "w");
    FILE *omegas = fopen("omegas.dat", "w");
    omega = 0.4;
    for(omega; omega<=1.0; omega+=0.2)
    {
        for(j=0; j<7; j++)
        {
            current = 100;
            E = omega*(1.4+j) - 0.1;
            while (current > acceptance)
            {
                E += 0.001/omega;
                if (j%2 == 0)
                {
                    psi[0] = 0;
                    psi_prime = 1;
                }
                else
                {
                    psi[0] = 1;
                    psi_prime = 0;
                }
                psi = solver();
                current = fabs(psi[N*5/10-1]);
            }
            printf("%f %f\n", omega, E);
            for(i = 0; i<N; i += 10)
                {
                    fprintf(functions, "%f\n", psi[i]);
                }
                fprintf(energies, "%f\n", E);
        }
        fprintf(omegas, "%f\n", omega);
    }
	return 0;
}

double *solver()
{
    int i = 0;
    double x = 0, U = 0;
    for(i = 1; i < N; i++)
    {
        U = 0.5*pow(omega*x, 2);
        x = x + dx;
        psi_prime = psi_prime + 2*(U-E)*psi[i-1]*dx;
        psi[i] = psi[i-1] + psi_prime*dx;
    }
    return psi;
}

