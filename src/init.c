#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "dg.h"
#include "dg1d.h"

void SetTestCaseData();
void InitCondShocktube(REAL x, REAL * U);
void InitCondBlast(REAL x, REAL * U);
void InitCondShuOsher(REAL x, REAL * U);

void SetTestCaseData()
{
   if(test_case == SOD)
   {
      printf("Setting SOD test case\n");
      xmin = 0.0;
      xmax = 1.0;
      XS   = 0.5;
      finaltime = 0.2;
      
      d_left = 1.0;   u_left = 0.0; p_left = 1.0;
      d_right= 0.125; u_right= 0.0; p_right= 0.1;
      bc_left = FREE;
      bc_right= FREE;
   }
   else if(test_case == BLAST)
   {
      printf("Setting BLAST test case\n");
      xmin = 0.0;
      xmax = 1.0;
      finaltime = 0.038;
      // reflecting bc
   }
   else if(test_case == LAX)
   {
      printf("Setting LAX test case\n");
      xmin = -5.0;
      xmax = +5.0;
      XS = 0.0;
      finaltime = 1.3;
      d_left  = 0.445;
      d_right = 0.5;
      
      u_left  = 0.698;
      u_right = 0.0;
      
      p_left  = 3.528;
      p_right = 0.571;
      bc_left = FREE;
      bc_right= FREE;
   }
   else if(test_case == LOWD)
   {
      printf("Setting LOWD test case\n");
      xmin = 0.0;
      xmax = 1.0;
      XS = 0.5;
      finaltime = 0.15;
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = -2.0;
      u_right =  2.0;
      
      p_left  = 0.4;
      p_right = 0.4;
      bc_left = FREE;
      bc_right= FREE;
   }
   else if(test_case == SHUOSHER)
   {
      printf("Setting SHUOSHER test case\n");
      xmin = -10.0;
      xmax =  5.0;
      finaltime = 1.8;
      bc_left = FREE;
      bc_right= FREE;
   }
   else
   {
      printf("Error in SetTestCaseData\n");
      exit(0);
   }
}

/* Allocate memory for main structure and set initial conditions */
CELL* Init()
{
   void ReadInput();
   void InitCondEuler(REAL, REAL *);
   UINT i, j, k, l;
   REAL dx, U[NVAR], v;
   CELL *cell;

   ReadInput();
   SetTestCaseData();

   // Note that PORD = degree + 1
   // Number of Gauss quadrature points
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
      cell[i].active = true;
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

      cell[i].U  = (REAL **) calloc(NVAR, sizeof(REAL *));
      cell[i].Re = (REAL **) calloc(NVAR, sizeof(REAL *));
      for(j = 0; j < NVAR; j++)
      {
         cell[i].U[j]  = (REAL *) calloc(cell[i].p, sizeof(REAL));
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
            cell[i].U[j][k] = 0.0;

      for(j = 0; j < cell[i].p; j++)
         for(k = 0; k < cell[i].ng; k++)
         {
            InitCondEuler(cell[i].xg[k], U);
            v = ShapeFun(cell[i].xg[k], &cell[i], j);
            for(l = 0; l < NVAR; l++)
               cell[i].U[l][j] += 0.5 * U[l] * v * wg[cell[i].ng - 1][k];
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
   if(fp == NULL)
   {
      printf("Error: Could not open inp.dat\n");
      exit(0);
   }
   
   fscanf(fp, "%s%lf", dummy, &cfl);
   fscanf(fp, "%s%d", dummy, &NC);
   fscanf(fp, "%s%d", dummy, &PORD);
   fscanf(fp, "%s%d", dummy, &NPLT);
   fscanf(fp, "%s%d", dummy, &FLUX);
   fscanf(fp, "%s%lf", dummy, &Mfact);
   fscanf(fp, "%s%d", dummy, &ALE);
   fscanf(fp, "%s%d", dummy, &test_case);
   fclose(fp);
   
   NF = NC + 1;
}

/* Initial condition for Euler equation */
void InitCondEuler(REAL x, REAL * U)
{
   REAL V[NVAR];
   
   if(test_case == SOD || test_case == LAX || test_case == LOWD)
      InitCondShocktube(x, V);
   else if(test_case == BLAST)
      InitCondBlast(x, V);
   else if(test_case == SHUOSHER)
      InitCondShuOsher(x, V);
   else
   {
      printf("Error: Unknown test case\n");
      exit(0);
   }
   
   U[0] = V[0];
   U[1] = V[0] * V[1];
   U[2] = V[2] / (GAMMA - 1.0) + 0.5 * V[0] * V[1] * V[1];
}

void InitCondShocktube(REAL x, REAL * V)
{
   if(x < XS)
   {
      V[0] = d_left;
      V[1] = u_left;
      V[2] = p_left;
   }
   else
   {
      V[0] = d_right;
      V[1] = u_right;
      V[2] = p_right;
   }
}

void InitCondBlast(REAL x, REAL * V)
{
   if(x < 0.1)
   {
      V[0] = 1.0;
      V[1] = 0.0;
      V[2] = 1000.0;
   }
   else if(x > 0.9)
   {
      V[0] = 1.0;
      V[1] = 0.0;
      V[2] = 100.0;
   }
   else
   {
      V[0] = 1.0;
      V[1] = 0.0;
      V[2] = 0.01;
   }
}

void InitCondShuOsher(REAL x, REAL * V)
{
   if(x < -4.0)
   {
      V[0] = 3.857143;
      V[1] = 2.629369;
      V[2] = 10.333333;
   }
   else
   {
      V[0] = 1 + 0.2*sin(5.0*x);
      V[1] = 0;
      V[2] = 1;
   }
}
