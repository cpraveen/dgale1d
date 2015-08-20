#include<stdio.h>
#include<stdlib.h>
#include "dg.h"
#include "dg1d.h"

int main()
{
   UINT iter;
   REAL time;
   CELL *cell;
   FACE *face;

   GaussInit();
   cell = Init();
   face = InitFaces(cell);
   
//   Result(cell); exit(0);

   cfl = cfl/(2*(PORD-1)+1);
   time = 0.0;
   iter = 0;

   printf("Beginning of iterations ...\n");
   while(time < finaltime)
   {
      // compute mesh velocity
      MeshVel(face);
      
      TimeStep(cell);
      if(time + dt > finaltime)
         dt = finaltime - time;
      
      // Assemble residual
      Flux (cell, face);
      // update solution
      Update (cell, face);
      // apply limiter
      Project (cell);

      // merge small cells, divide large
      
      time += dt;
      ++iter;
      printf("%8d  %16.6e %16.6e %16.6e %16.6e\n", iter, dt, time, dxmin, dxmax);
   }
   Result(cell);

   return 0;
}
