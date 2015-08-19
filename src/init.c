#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "dg.h"
#include "dg1d.h"

/* Allocate memory for main structure and set initial conditions */
CELL* Init()
{
   void ReadInput();
   void InitCondEuler(REAL, REAL *);
   UINT i, j, k, l;
   REAL dx, U[NVAR], v;
   CELL *cell;

   ReadInput();

   /* Coefficients for RK3 */
   ark[0] = 0.0;
   ark[1] = 3.0 / 4.0;
   ark[2] = 1.0 / 3.0;

   brk[0] = 1.0;
   brk[1] = 1.0 / 4.0;
   brk[2] = 2.0 / 3.0;

   // PORD = degree + 1
   NG = PORD;

   printf("Allocating memory and setting initial condition ...\n");

   cell = (CELL *) calloc(NC, sizeof(CELL));
   if(cell == NULL)
   {
      printf("Init: Could not allocate cell\n");
      exit(0);
   }
   
   dx = (xmax - xmin) / NC;
   printf("No of cells = %d\n", NC);
   printf("No of faces = %d\n", NF);
   printf("         dx = %f\n", dx);

   // Initialize cells
   for(i = 0; i < NC; i++)
   {
      cell[i].xl = xmin + i * dx;
      cell[i].xr = cell[i].xl + dx;
      cell[i].x = 0.5 * (cell[i].xl + cell[i].xr);
      cell[i].h = cell[i].xr - cell[i].xl;
      
      cell[i].wl = 0.0;
      cell[i].wr = 0.0;

      cell[i].p = PORD;

      cell[i].ng = NG;
      cell[i].xg = (REAL *) calloc(cell[i].ng, sizeof(REAL));
      GaussPoints(&cell[i]);

      cell[i].Un = (REAL **) calloc(NVAR, sizeof(REAL *));
      cell[i].Uo = (REAL **) calloc(NVAR, sizeof(REAL *));
      cell[i].Re = (REAL **) calloc(NVAR, sizeof(REAL *));
      for(j = 0; j < NVAR; j++)
      {
         cell[i].Un[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
         cell[i].Uo[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
         cell[i].Re[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
      }
      
      // Needs changes for periodic bc
      if(i == 0)
      {
         cell[i].lcell = &cell[i];
         cell[i].rcell = &cell[i+1];
      }
      else if(i == NC-1)
      {
         cell[i].lcell = &cell[i-1];
         cell[i].rcell = &cell[i];
      }
      else
      {
         cell[i].lcell = &cell[i-1];
         cell[i].rcell = &cell[i+1];
      }
   }
   
   /* Set initial condition by L2 projection */
   for(i = 0; i < NC; i++)
   {

      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++)
            cell[i].Un[j][k] = 0.0;

      for(j = 0; j < cell[i].p; j++)
         for(k = 0; k < cell[i].ng; k++)
         {
            InitCondEuler(cell[i].xg[k], U);
            v = ShapeFun(cell[i].xg[k], &cell[i], j);
            for(l = 0; l < NVAR; l++)
               cell[i].Un[l][j] += 0.5 * U[l] * v * wg[cell[i].ng - 1][k];
         }
   }

   return cell;
}

FACE* InitFaces(CELL *cell)
{
   UINT i;
   REAL dx;
   FACE *face;
   
   face = (FACE *) calloc(NF, sizeof(FACE));
   if(face == NULL)
   {
      printf("Init: Could not allocate face\n");
      exit(0);
   }
   
   dx = (xmax - xmin) / NC;

   // Initialize faces
   face[0].x  = xmin;
   face[0].w  = 0.0;
   face[0].lcell = &cell[0]; // change this for periodic bc
   face[0].rcell = &cell[0];
   face[0].active = true;
   for(i = 1; i < NF-1; i++)
   {
      face[i].x = xmin + i * dx;
      face[i].w = 0.0;
      face[i].lcell = &cell[i-1];
      face[i].rcell = &cell[i];
      face[i].active = true;
   }
   face[NF-1].x = xmax;
   face[NF-1].w = 0.0;
   face[NF-1].lcell = &cell[NC-1]; // change this for periodic bc
   face[NF-1].rcell = &cell[NC-1];
   face[NF-1].active = true;

   return face;
}

void ReadInput()
{
   FILE *fp;
   char dummy[100];
   fp = fopen("inp.dat", "r");
   if(fp == NULL) {
      printf("Error: Could not open inp.dat\n");
      exit(0);
   }
   
   fscanf(fp, "%s%lf", dummy, &cfl);
   fscanf(fp, "%s%lf", dummy, &finaltime);
   fscanf(fp, "%s%d", dummy, &NC);
   fscanf(fp, "%s%d", dummy, &PORD);
   fscanf(fp, "%s%d", dummy, &NPLT);
   fscanf(fp, "%s%d", dummy, &FLUX);
   fscanf(fp, "%s%lf", dummy, &Mfact);
   fscanf(fp, "%s%lf%lf", dummy, &xmin, &xmax);
   fscanf(fp, "%s%lf", dummy, &XS);
   fscanf(fp, "%s%lf%lf%lf", dummy, &d_left, &u_left, &p_left);
   fscanf(fp, "%s%lf%lf%lf", dummy, &d_right, &u_right, &p_right);
   fclose(fp);
   
   NF = NC + 1;
}

/* Initial condition for Burgers equation */
REAL InitCondBurger(REAL x)
{
   if(x < 0.5)
      return 1.0;
   else
      return 0.0;
}

/* Initial condition for Euler equation */
void InitCondEuler(REAL x, REAL * U)
{
   REAL d, u, p;

   if(x < XS) {
      d = d_left;
      u = u_left;
      p = p_left;
   }
   else {
      d = d_right;
      u = u_right;
      p = p_right;
   }

   U[0] = d;
   //U[0] = sin(M_PI*x);
   U[1] = d * u;
   U[2] = p / (GAMMA - 1.0) + 0.5 * d * u * u;
}
