#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 100000;
double dx = 0.0001, E, omega;
double *psi;
double psi_prime;

double *solver();
double *leapfrog();
double *kutta();
double force(double x, double y);

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

double *kutta()
{
    double k1, k2, k3, k4, k, l1, l2, l3, l4, l, x=0;
    int i = 0;
    double temp;
    temp = psi_prime;
    for(i = 1; i < N; i++)
    {
        l1 = force(x, psi[i-1]);
        l2 = force(x + 0.5*dx, psi[i-1] + 0.5*l1);
        l3 = force(x + 0.5*dx, psi[i-1] + 0.5*l2);
        l4 = force(x + dx, psi[i-1] + l3);
        l = (l1 + 2*l2 + 2*l3 + l4)*dx/6.0;
        psi_prime += l;
        
        k1 = psi_prime;
        k2 = k1 + 0.5*dx*k1;
        k3 = k1 + 0.5*dx*k2;
        k4 = k1 + dx*k3;
        k = (k1 + 2*k2 + 2*k3 + k4)*dx/6.0;
       
        psi[i] = psi[i-1] + k;
        x += dx;
    }
    return psi;
}
double *leapfrog()
{
    double half, x = 0;
    int i = 0, U = 0;
    for(i = 1; i<N; i++)
    {
        half = psi_prime + 0.5*force(x, psi[i-1])*dx;
        psi[i] = psi[i-1] + dx*half;
        psi_prime = half + 0.5*force(x, psi[i])*dx;
        x += dx;
    }
    return psi;
}

double force(double x, double y)
{
    double U;
    U = 0.5*pow(omega*x, 2);
    return 2*(U - E)*y;
}

//double *adams()
//{
    //int i = 0;
    //double x = 0, U = 0;
    //for(i = 2; i<N; i++)
    //{
        //psi[i] = psi[i-1] + 1.5*dx*force(x, psi[i-1]) - 0.5*dx*force(x-dx, psi[i-2])
        //x +dx 
    //}
    //return psi;
//}

