#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

void ApplyPositivityLimiter(CELL *cell);

void Project(CELL * cell)
{
   REAL minmod(REAL, REAL, REAL);
   REAL minmod2(REAL, REAL, REAL, REAL, REAL);
   UINT i, j, k;
   REAL u[NVAR], ux[NVAR], uxb[NVAR], dul[NVAR], dur[NVAR], R[NVAR][NVAR],
      Ri[NVAR][NVAR], fact;

   fact = sqrt(3.0);

   for(i = 0; i < NC; i++)
      if(cell[i].active)
      {
         for(j = 0; j < NVAR; j++)
         {
            dul[j] = cell[i].U[j][0] - cell[i].lcell->U[j][0];
            dur[j] = cell[i].rcell->U[j][0] - cell[i].U[j][0];
            u[j] = cell[i].U[j][0];
            ux[j] = fact * cell[i].U[j][1];
         }
         
         EigMat(u, R, Ri);
         Multi(Ri, ux);
         Multi(Ri, dul);
         Multi(Ri, dur);
         for(j = 0; j < NVAR; j++)
         {
            if(fabs(ux[j]) <= Mfact * cell[i].h * cell[i].h)
               uxb[j] = ux[j];
            else
               uxb[j] = minmod(ux[j], dul[j], dur[j]);
         }
         Multi(R, uxb);
         
         for(j = 0; j < NVAR; j++)
         {
            uxb[j] = uxb[j] / fact;
            if(fabs(cell[i].U[j][1] - uxb[j]) > 1.0e-6)
            {
               cell[i].U[j][1] = uxb[j];
               for(k = 2; k < cell[i].p; k++)
                  cell[i].U[j][k] = 0.0;
            }
            
         }
         if(pos_lim == TRUE) ApplyPositivityLimiter(&cell[i]);
      }
}

/* minmod limiter function */
REAL minmod(REAL a, REAL b, REAL c)
{
   REAL sgn, m;

   if(a * b <= 0.0 || b * c <= 0.0)
      return 0.0;

   sgn = (a > 0.0) ? 1.0 : -1.0;
   a = fabs(a);
   b = fabs(b);
   c = fabs(c);
   m = (a < b) ? a : b;
   m = (c < m) ? c : m;
   return sgn * m;

}

/* Eigenvector matrix */
void EigMat(REAL * U, REAL R[][3], REAL Ri[][3])
{
   REAL d, v, p, c, h, f, g1, g2;

   g1 = GAMMA - 1.0;
   g2 = g1 / 2.0;

   d = U[0];
   v = U[1] / d;
   p = (GAMMA - 1.0) * (U[2] - 0.5 * d * v * v);
   c = sqrt(GAMMA * p / d);
   h = c * c / g1 + 0.5 * v * v;
   f = d / c / 2.0;

   /* Inverse eigenvector-matrix */
   Ri[0][0] = 1.0 - g2 * v * v / c / c;
   Ri[1][0] = (g2 * v * v - v * c) / d / c;
   Ri[2][0] = -(g2 * v * v + v * c) / d / c;

   Ri[0][1] = g1 * v / c / c;
   Ri[1][1] = (c - g1 * v) / d / c;
   Ri[2][1] = (c + g1 * v) / d / c;

   Ri[0][2] = -g1 / c / c;
   Ri[1][2] = g1 / d / c;
   Ri[2][2] = -g1 / d / c;

   /* Eigenvector matrix */
   R[0][0] = 1.0;
   R[1][0] = v;
   R[2][0] = v * v / 2.0;

   R[0][1] = f;
   R[1][1] = (v + c) * f;
   R[2][1] = (h + v * c) * f;

   R[0][2] = -f;
   R[1][2] = -(v - c) * f;
   R[2][2] = -(h - v * c) * f;

}

/* Multiply matrix R and vector U */
void Multi(REAL R[][3], REAL * U)
{
   UINT i, j;
   REAL Ut[NVAR];

   for(i = 0; i < NVAR; i++)
      Ut[i] = U[i];

   for(i = 0; i < NVAR; i++)
   {
      U[i] = 0.0;
      for(j = 0; j < NVAR; j++)
         U[i] += R[i][j] * Ut[j];
   }
}

// Apply positivity limiter in cell
void ApplyPositivityLimiter(CELL *cell)
{
   UINT i, j;
   const REAL eps = 1.0e-13;
   REAL theta1, theta2, rho_min, rat, pre, t1, t2, t, a1, b1, c1, D;
   REAL drho, dm, dE;
   REAL **U;
   
   // First order scheme, nothing to do
   if(cell->p == 1) return;
   
   U  = (REAL **) calloc(cell->ngll, sizeof(REAL *));
   for(i = 0; i < cell->ngll; i++)
      U[i]  = (REAL *) calloc(NVAR, sizeof(REAL));

   // First, limit density
   
   // Compute solution at GLL nodes
   UatGLL(cell, U);

   // Find minimum value of density
   rho_min = 1.0e20;
   for(i=0; i<cell->ngll; ++i)
      rho_min = MIN(rho_min, U[i][0]);
   
   rat = fabs(cell->U[0][0] - eps)/(fabs(cell->U[0][0] - rho_min) + 1.0e-14);
   theta1 = MIN(1.0, rat);
   
   // Apply limiter, dont change mean value
   for(i=0; i<NVAR; ++i)
      for(j=1; j<cell->p; ++j)
         cell->U[i][j] *= theta1;
   
   // Now limit the pressure

   // Compute solution at GLL nodes
   UatGLL(cell, U);
   
   theta2 = 1.0;
   for(i=0; i < cell->ngll; ++i)
   {
      pre = (GAMMA-1.0) * (U[i][2] - 0.5 * U[i][1] * U[i][1] / U[i][0]);
      if(pre < eps)
      {
         drho = U[i][0] - cell->U[0][0];
         dm   = U[i][1] - cell->U[1][0];
         dE   = U[i][2] - cell->U[2][0];
         a1 = 2.0 * drho * dE - dm * dm;
         b1 = 2.0 * drho * (cell->U[2][0] - eps/(GAMMA-1.0))
              + 2.0 * cell->U[0][0] * dE
              - 2.0 * cell->U[1][0] * dm;
         c1 = 2.0 * cell->U[0][0] * cell->U[2][0]
              - pow(cell->U[1][0],2)
              - 2.0 * eps * cell->U[0][0]/(GAMMA-1.0);
         // Divide by a1 to avoid round-off error
         b1 /= a1; c1 /= a1;
         D = sqrt( fabs(b1*b1 - 4.0*c1) );
         t1 = 0.5*(-b1 - D);
         t2 = 0.5*(-b1 + D);
         if(t1 > -1.0e-12 && t1 < 1.0 + 1.0e-12)
            t = t1;
         else if(t2 > -1.0e-12 && t2 < 1.0 + 1.0e-12)
            t = t2;
         else
         {
            printf("Fatal error\n");
            printf("Mean rho = %e", cell->U[0][0]);
            printf("pre at gll point = %e\n", pre);
            printf("t1 = %e, t2 = %e\n", t1, t2);
            exit(0);
         }
         t = MIN(1.0, t);
         t = MAX(0.0, t);
         // Need t < 1.0. If t==1 upto machine precision
         // then we are suffering from round off error.
         // In this case we take the cell average value, t=0.
         if(fabs(1.0-t) < 1.0e-14) t = 0.0;
         theta2 = MIN(theta2, t);
      }
   }

   // Apply limiter, dont change mean value
   for(i=0; i<NVAR; ++i)
      for(j=1; j<cell->p; ++j)
         cell->U[i][j] *= theta2;
   
   // Release memory
   for(i = 0; i < cell->ngll; i++)
      free(U[i]);
   free(U);
}
