#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int N = 10000, Z = 2, l=0;
double dx, *V, *U;

double *linspace(double min, double max, int N);
void initial_potential();
double integrate(double *function, double dx, int N);
double *cal_probability(double *function);
void solver(double *function, double *derivative, double epsilon, int l);
double seaker(double *function, double *derivative);
double *cal_charge(double *probability);
void cal_potential(double V_0, double *charge);

int main(int argc, char **argv)
{
    int i;
    double der, epsilon;
    double *R = malloc(N*sizeof(double));
    double *R_prime = malloc(N*sizeof(double));
    double *P, *charge;
    V = malloc(N*sizeof(double));
	U = linspace(0.001, 30, N);
    initial_potential();

    R[0] = 1;
    R_prime[0] = -0.99;
    epsilon = seaker(R, R_prime);
    P = cal_probability(R);
    charge = cal_charge(P);
    cal_potential(0, charge);
    
    FILE *output = fopen("data.dat", "w");

    for(i=0; i<N; i++)
    {
        fprintf(output, "%f %f %f %f %f\n", U[i], R[i], P[i], charge[i], V[i]);
    }
    fclose(output);
	return 0;
}

void initial_potential()
{
    int i;
    double coeff, exponent;
    for(i = 0; i<N; i++)
    {
        coeff = -2/(Z*U[i]);
        exponent = exp(-U[i]);
        V[i] = coeff*(1+(Z-1)*exponent);
    }
}

double *cal_charge(double *probability)
{
    int i;
    double *charge = malloc(N*sizeof(double));
    for(i=0; i<N; i++)
    {
        charge[i] = Z-(Z-1)*integrate(probability, dx, i+1);
    }
    return charge;
}

void cal_potential(double V_0, double *charge)
{
    int i;
    double *delta = malloc(N*sizeof(double));
    for(i=0; i<N; i++)
    {
        delta[i] = -(2/Z)*(charge[i]/(U[i]*U[i]));
    }
    V[0] = delta[0] + V_0;
    for(i=1; i<N; i++)
    {
        V[i] = integrate(delta, dx, i+1);
    }
    free(delta);
}
        
double *cal_probability(double *function)
{
    int i;
    double *P = malloc(N*sizeof(double)), norm;
    for(i=0; i<N; i++)
    {
        P[i] = pow(U[i]*function[i], 2);
    }
    
    norm = integrate(P, dx, N);
    for(i=0; i<N; i++)
    {
        P[i] *= 1/norm;
    }
    return P;
}

double integrate(double *function, double dx, int N)
{
    int i;
    double integral = 0;
    for(i=0; i<N; i++)
    {
        integral += function[i];
    }
    return integral*dx;
}
    
double seaker(double *function, double *derivative)
{
    double energy, last, de, current;
    
    energy = -0.5;
    de = 0.01;
    last = 0;
    current = 1;
    while((fabs(current) >= 1E-3) && (de > 1E-9))
    {
        energy += de;
        solver(function, derivative, energy, l);
        current = function[N-1];
        if (current*last < 0)
        {
            energy += -de;
            de *= 0.1;
        }
        last = current;
    }
    return energy;
}

double *linspace(double min, double max, int N)
{
    int i;
    dx = (max-min)/(N-1);
    double *x = malloc(N*sizeof(double));
    
    x[0] = min;
    for(i=0; i<(N-1); i++)
    {
        x[i+1] = x[i] + dx;
    }
    return x;
}

void solver(double *function, double *derivative, double epsilon, int l)
{
    int L = l*(l+1), i = 0;
    double der;
    
    for(i=0; i<(N-1); i++)
    {
        der = -2*derivative[i]/U[i] - function[i]*(epsilon - V[i] - L/(U[i]*U[i]));
        derivative[i+1] = derivative[i] + der*dx;
        function[i+1] = function[i] + derivative[i+1]*dx;
    }
}
