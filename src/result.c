#include <stdio.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

void Result(CELL * cell)
{
   FILE *fp1, *fp2, *fp3;
   UINT i, j;
   REAL dx, x, U[NVAR], d, u, p, c, m, w;
   
   fp1 = fopen("avg", "w");
   fp2 = fopen("sol", "w");
   fp3 = fopen("h", "w");

   for(i = 0; i < NC; i++)
      if(cell[i].active)
      {
         // average solution
         d = cell[i].U[0][0];
         u = cell[i].U[1][0] / d;
         p = (GAMMA - 1.0) * (cell[i].U[2][0] - 0.5 * d * u * u);
         c = sqrt(GAMMA * p / d);
         m = u / c;
         fprintf(fp1, "%e %e %e %e %e\n", cell[i].x, d, u, p, m);

         // more detailed solution evaluated at NPLT points inside cell
         if(NPLT == 1)
            dx = 0.5 * cell[i].h;
         else
            dx = cell[i].h / (NPLT - 1);
         for(j = 0; j < NPLT; j++)
         {
            if(NPLT == 1)
               x = cell[i].xl + dx;
            else
               x = cell[i].xl + dx * j;
            Uvect(&cell[i], x, U);
            d = U[0];
            u = U[1] / d;
            p = (GAMMA - 1.0) * (U[2] - 0.5 * d * u * u);
            c = sqrt(GAMMA * p / d);
            m = u / c;
            w = GetMeshVel(&cell[i], x);
            fprintf(fp2, "%f %f %f %f %f %f\n", x, d, u, p, m, w);
         }
         fprintf(fp2, "\n");
         
         // mesh size
         fprintf(fp3, "%d %e %e\n", i, cell[i].x, cell[i].h);
      }
   fclose(fp1);
   fclose(fp2);
   fclose(fp3);
}
