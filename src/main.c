#include<stdio.h>
#include<stdlib.h>
#include "dg.h"
#include "dg1d.h"

int main()
{
   CELL* Init();
   FACE* InitFaces(CELL*);
   void TimeStep(CELL *);
   void SaveSol(CELL *);
   void Flux(CELL *, FACE * face);
   void Update(UINT, CELL *);
   void Project(CELL *);
   void Result(CELL *);

   UINT iter, rk;
   REAL time;
   CELL *cell;
   FACE *face;

   NVAR = 3;                    /* Number of variables */
   RK = 3;                      /* Number of Runge-Kutta stages */

   GaussInit();
   cell = Init();
   face = InitFaces(cell);
   printf("aaa\n");
   printf("%e\n", cell[0].x);
   printf("bbb\n");
   
//   Result(cell); exit(0);

   cfl = cfl / PORD;
   time = 0.0;
   iter = 0;

   printf("Beginning of iterations ...\n");
   while(time < finaltime)
   {
      SaveSol(cell);
      TimeStep(cell);
      if(time + dt > finaltime)
         dt = finaltime - time;
      for(rk = 0; rk < RK; rk++)
      {
         Flux(cell, face);
         Update(rk, cell);
         Project(cell);
      }
      time += dt;
      ++iter;
      printf("%8d  %18.6e %18.6e\n", iter, dt, time);
   }
   Result(cell);

   return 0;
}
