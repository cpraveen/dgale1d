#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

REAL RoeEigVal(REAL *, REAL *);
REAL MaxEigVal(REAL *, REAL *, REAL);

//-----------------------------------------------------------------------------
/* Compute flux for 1-d euler equation given conserved vector */
//-----------------------------------------------------------------------------
void EulerFlux(REAL * U, REAL w, REAL * flux)
{
   REAL p;

   p = (GAMMA - 1.0) * (U[2] - 0.5 * U[1] * U[1] / U[0]);

   flux[0] = U[1] - w*U[0];
   flux[1] = p + U[1] * U[1] / U[0] - w*U[1];
   flux[2] = (U[2] + p) * U[1] / U[0] - w*U[2];
}
//-----------------------------------------------------------------------------
// Compute flux Jacobian
//-----------------------------------------------------------------------------
void Jacobian(REAL *U, REAL A[][3])
{
   REAL v  = U[1]/U[0];
   REAL v2 = v * v;
   REAL e  = U[2]/U[0];
   
   A[0][0] = 0;                           A[0][1] = 1.0;                        A[0][2] = 0;
   A[1][0] = 0.5*(GAMMA-3)*v2;            A[1][1] = (3-GAMMA)*v;                A[1][2] = GAMMA-1;
   A[2][0] = (-GAMMA*e + (GAMMA-1)*v2)*v; A[2][1] = GAMMA*e - 1.5*(GAMMA-1)*v2; A[2][2] = GAMMA*v;
}
//-----------------------------------------------------------------------------
/* Roe flux for 1-d euler equations */
//-----------------------------------------------------------------------------
void RoeFlux (REAL * Ul, REAL * Ur, REAL w, REAL * flux)
{
   REAL left[NVAR], right[NVAR];
   UINT i, j;

   left[0] = Ul[0];
   left[1] = Ul[1] / Ul[0];
   left[2] = (GAMMA-1) * ( Ul[2] - 0.5 * left[0] * left[1] * left[1]);

   right[0] = Ur[0];
   right[1] = Ur[1] / Ur[0];
   right[2] = (GAMMA-1) * ( Ur[2] - 0.5 * right[0] * right[1] * right[1]);

   REAL fl  = sqrt(left[0]);
   REAL fr  = sqrt(right[0]);
   REAL u   = (fl*left[1] + fr*right[1])/(fl + fr);

   REAL Hl  = GAMMA*left[2]/left[0]/(GAMMA-1.0) + 0.5*pow(left[1],2);
   REAL Hr  = GAMMA*right[2]/right[0]/(GAMMA-1.0) + 0.5*pow(right[1],2);

   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] - w*Ul[0] + right[0]*right[1] - w*Ur[0]);
   flux[1] = 0.5*(left[2]  + left[0]  * pow(left[1],2)  - w*Ul[1] +
                  right[2] + right[0] * pow(right[1],2) - w*Ur[1]);
   flux[2] = 0.5*(Hl*left[0]*left[1] - w*Ul[2] + Hr*right[0]*right[1] - w*Ur[2]);

   // Add conservative dissipation
   REAL H = (fl*Hl + fr*Hr)/(fl + fr);
   REAL a = sqrt((GAMMA-1.0)*(H - 0.5*u*u));
   REAL R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;

   REAL Lambda[] = { fabs(u-w-a), fabs(u-w), fabs(u-w+a)};

   REAL dU[] = {
      Ur[0] - Ul[0],
      Ur[1] - Ul[1],
      Ur[2] - Ul[2]
   };

   REAL aa[3];
   aa[1] = (GAMMA-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];

   for(i=0; i<3; ++i)
      for(j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
}

//-----------------------------------------------------------------------------
/* Lax-Friedrichs flux for 1-d euler equations */
//-----------------------------------------------------------------------------
void LFFlux(REAL * Ul, REAL * Ur, REAL w, REAL * flux)
{
   UINT i;
   REAL Fl[NVAR], Fr[NVAR], lam;

   EulerFlux(Ul, w, Fl);
   EulerFlux(Ur, w, Fr);
   //lam = RoeEigVal (Ul, Ur);
   lam = MaxEigVal(Ul, Ur, w);

   for(i = 0; i < NVAR; i++)
      flux[i] = 0.5 * (Fl[i] + Fr[i] - lam * (Ur[i] - Ul[i]));
}

//-----------------------------------------------------------------------------
/* Maximum eigenvalue of left and right states */
//-----------------------------------------------------------------------------
REAL MaxEigVal(REAL * Ul, REAL * Ur, REAL w)
{
   REAL dl, ul, pl, al, ll, dr, ur, pr, ar, lr;

   dl = Ul[0];
   ul = Ul[1] / dl;
   pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
   al = sqrt(GAMMA * pl / dl);
   ll = fabs(ul-w) + al;

   dr = Ur[0];
   ur = Ur[1] / dr;
   pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
   ar = sqrt(GAMMA * pr / dr);
   lr = fabs(ur-w) + ar;

   if(ll > lr)
      return ll;
   else
      return lr;
}

//-----------------------------------------------------------------------------
/* Returns maximum eigenvalue abs(u) + a based on Roe-average */
//-----------------------------------------------------------------------------
REAL RoeEigVal(REAL * Ul, REAL * Ur)
{
   REAL dl, ul, pl, hl, dr, ur, pr, hr, d1, d2, d3, u, h, a;

   dl = Ul[0];
   ul = Ul[1] / dl;
   pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
   hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;

   dr = Ur[0];
   ur = Ur[1] / dr;
   pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
   hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;

   d1 = sqrt(dl);
   d2 = sqrt(dr);
   d3 = 1.0 / (d1 + d2);

   u = (d1 * ul + d2 * ur) * d3;
   h = (d1 * hl + d2 * hr) * d3;
   a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u));

   return fabs(u) + a;
}

//-----------------------------------------------------------------------------
/* Characteristic Lax-Friedrichs flux for 1-d euler equations */
//-----------------------------------------------------------------------------
//void LFCFlux(REAL * Ul, REAL * Ur, REAL * flux)
//{
//   void RoeAverage(REAL *, REAL *, REAL *);
//   UINT i;
//   REAL U[NVAR], Fl[NVAR], Fr[NVAR], R[NVAR][NVAR], Ri[NVAR][NVAR], lam[NVAR];
//   REAL d, u, p, a;
//
//   EulerFlux(Ul, Fl);
//   EulerFlux(Ur, Fr);
//   RoeAverage(Ul, Ur, U);
//   EigMat(U, R, Ri);
//   Multi(Ri, Fl);
//   Multi(Ri, Fr);
//   Multi(Ri, Ul);
//   Multi(Ri, Ur);
//
//   /* Eigenvalues based on Roe-averaged state */
//   d = U[0];
//   u = U[1] / d;
//   p = (GAMMA - 1.0) * (U[2] - 0.5 * d * u * u);
//   a = sqrt(GAMMA * p / d);
//   lam[0] = fabs(u);
//   lam[1] = fabs(u + a);
//   lam[2] = fabs(u - a);
//
//   for(i = 0; i < NVAR; i++)
//      flux[i] = 0.5 * (Fl[i] + Fr[i] - lam[i] * (Ur[i] - Ul[i]));
//
//   Multi(R, flux);
//}

//-----------------------------------------------------------------------------
/* Returns maximum eigenvalue abs(u) + a based on Roe-average */
//-----------------------------------------------------------------------------
void RoeAverage(REAL * Ul, REAL * Ur, REAL * U)
{
   REAL dl, ul, pl, hl, dr, ur, pr, hr, d1, d2, d3, d, u, h, p;

   dl = Ul[0];
   ul = Ul[1] / dl;
   pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
   hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;

   dr = Ur[0];
   ur = Ur[1] / dr;
   pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
   hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;

   d1 = sqrt(dl);
   d2 = sqrt(dr);
   d3 = 1.0 / (d1 + d2);

   u = (d1 * ul + d2 * ur) * d3;
   h = (d1 * hl + d2 * hr) * d3;
   d = sqrt(dl * dr);
   p = (GAMMA - 1.0) * (h - 0.5 * u * u) * d / GAMMA;

   U[0] = d;
   U[1] = d * u;
   U[2] = p / (GAMMA - 1.0) + 0.5 * d * u * u;

}

//-----------------------------------------------------------------------------
/* ECUSP flux of Jameson */
//-----------------------------------------------------------------------------
//void ECUSPFlux(REAL * Ul, REAL * Ur, REAL * flux)
//{
//   UINT i;
//   REAL dl, ul, pl, hl, dr, ur, pr, hr, d1, d2, d3, u, h, a;
//   REAL m, m1, alpha, beta, du, Fp[3], Ua[3], Fl[3], Fr[3], dv;
//   REAL EPS = 0.01;
//
//   EulerFlux(Ul, Fl);
//   EulerFlux(Ur, Fr);
//
//   dl = Ul[0];
//   ul = Ul[1] / dl;
//   pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
//   hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;
//
//   dr = Ur[0];
//   ur = Ur[1] / dr;
//   pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
//   hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;
//
//   /* Roe average values */
//   d1 = sqrt(dl);
//   d2 = sqrt(dr);
//   d3 = 1.0 / (d1 + d2);
//
//   u = (d1 * ul + d2 * ur) * d3;
//   h = (d1 * hl + d2 * hr) * d3;
//   a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u));
//   m = u / a;
//
//   alpha = fabs(m);
//
//   if(alpha >= 1.0)
//      beta = (m > 0.0) ? 1.0 : -1.0;
//   else if(m >= 0.0) {
//      m1 = 2.0 * m - 1.0;
//      beta = (m1 > 0.0) ? m1 : 0.0;
//   }
//   else {
//      m1 = 2.0 * m + 1.0;
//      beta = (m1 < 0.0) ? m1 : 0.0;
//   }
//
//   if(alpha < EPS)
//      alpha = 0.5 * (EPS + m * m / EPS);
//
//   du = ur - ul;
//   for(i = 0; i < 3; i++)
//      Ua[i] = 0.5 * (Ul[i] + Ur[i]);
//
//   Fp[0] = 0.0;
//   Fp[1] = pr - pl;
//   Fp[2] = pr * ur - pl * ul;
//
//   for(i = 0; i < 3; i++) {
//      dv = alpha * a * (Ur[i] - Ul[i]) + beta * du * Ua[i] + beta * Fp[i];
//      flux[i] = 0.5 * (Fl[i] + Fr[i] - dv);
//   }
//
//}

//-----------------------------------------------------------------------------
/* HLLC flux */
//-----------------------------------------------------------------------------
//void HLLCFlux(REAL * Ul, REAL * Ur, REAL * flux)
//{
//   UINT i;
//   REAL dl, ul, pl, hl, dr, ur, pr, hr;
//   REAL Fl[3], Fr[3];
//   REAL al, ar, da, aa, ps, us, ql, qr, Sl, Sr, f, Us[3];
//
//   dl = Ul[0];
//   ul = Ul[1] / dl;
//   pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
//   al = sqrt(GAMMA * pl / dl);
//   hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;
//
//   dr = Ur[0];
//   ur = Ur[1] / dr;
//   pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
//   ar = sqrt(GAMMA * pr / dr);
//   hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;
//
//   da = 0.5 * (dl + dr);
//   aa = 0.5 * (al + ar);
//   ps = 0.5 * (pl + pr - (ur - ul) * da * aa);
//   us = 0.5 * (ul + ur - (pr - pl) / da / aa);
//
//   if(ps <= pl)
//      ql = 1.0;
//   else
//      ql = sqrt(1.0 + 0.5 * (GAMMA + 1.0) * (ps / pl - 1.0) / GAMMA);
//
//   if(ps <= pr)
//      qr = 1.0;
//   else
//      qr = sqrt(1.0 + 0.5 * (GAMMA + 1.0) * (ps / pr - 1.0) / GAMMA);
//
//   Sl = ul - al * ql;
//   Sr = ur + ar * qr;
//
//   if(0.0 <= Sl)
//      EulerFlux(Ul, flux);
//   else if(Sl <= 0.0 && 0.0 <= us) {
//      EulerFlux(Ul, Fl);
//      f = dl * (Sl - ul) / (Sl - us);
//      Us[0] = 1.0;
//      Us[1] = us;
//      Us[2] = Ul[2] / dl + (us - ul) * (us + pl / dl / (Sl - ul));
//      for(i = 0; i < 3; i++)
//         flux[i] = Fl[i] + Sl * f * (Us[i] - Ul[i]);
//   }
//   else if(us <= 0.0 && 0.0 <= Sr) {
//      EulerFlux(Ur, Fr);
//      f = dr * (Sr - ur) / (Sr - us);
//      Us[0] = 1.0;
//      Us[1] = us;
//      Us[2] = Ur[2] / dr + (us - ur) * (us + pr / dr / (Sr - ur));
//      for(i = 0; i < 3; i++)
//         flux[i] = Fr[i] + Sr * f * (Us[i] - Ur[i]);
//   }
//   else
//      EulerFlux(Ur, flux);
//
//
//}

//-----------------------------------------------------------------------------
/* AUSMDV flux */
//-----------------------------------------------------------------------------
//void AUSMDVFlux(REAL * ql, REAL * qr, REAL * fl)
//{
//
//   REAL rho_left, rho_right, vel_left, vel_right,
//      pre_left, pre_right, c_left, c_right, c_mean,
//      h_left, h_right, amach_l, amach_r,
//      alpha_l, alpha_r, temp1, temp2, u_lp, u_rm,
//      p_lp, p_rm, ak_half, p_half, s_half, flux_v1,
//      flux_d1, flux_p, flux_n, p_min, a_min, gam, gam1;
//
//   gam = 1.4;
//   gam1 = gam - 1.0;
//
//   rho_left = ql[0];
//   vel_left = ql[1] / ql[0];
//   pre_left = gam1 * (ql[2] - 0.5 * ql[1] * ql[1] / ql[0]);
//
//   rho_right = qr[0];
//   vel_right = qr[1] / qr[0];
//   pre_right = gam1 * (qr[2] - 0.5 * qr[1] * qr[1] / qr[0]);
//
//   c_left = sqrt(gam * pre_left / rho_left);
//   h_left = c_left * c_left / gam1 + 0.5 * vel_left * vel_left;
//
////......right state 
//   c_right = sqrt(gam * pre_right / rho_right);
//   h_right = c_right * c_right / gam1 + 0.5 * vel_right * vel_right;
//
////......mean speed of sound
//   c_mean = (c_left > c_right) ? c_left : c_right;
//
////......left and right mach no.
//   amach_l = fabs(vel_left) / c_mean;
//   amach_r = fabs(vel_right) / c_mean;
//
////......alpha
//   temp1 = pre_left / rho_left;
//   temp2 = pre_right / rho_right;
//
//   alpha_l = 2 * temp1 / (temp1 + temp2);
//   alpha_r = 2 * temp2 / (temp1 + temp2);
//
////......vel_left+
//   temp1 = vel_left + c_mean;
//   temp2 = vel_left + fabs(vel_left);
//   if(amach_l <= 1.0)
//      u_lp = alpha_l * (0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
//         + 0.5 * temp2;
//   else
//      u_lp = 0.5 * temp2;
//
////......vel_right-
//   temp1 = vel_right - c_mean;
//   temp2 = vel_right - fabs(vel_right);
//   if(amach_r <= 1.0)
//      u_rm = alpha_r * (-0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
//         + 0.5 * temp2;
//   else
//      u_rm = 0.5 * temp2;
//
////......pre_left+
//   temp1 = vel_left / c_mean + 1.0;
//   temp2 = 2.0 - vel_left / c_mean;
//   if(amach_l <= 1.0)
//      p_lp = 0.25 * pre_left * temp1 * temp1 * temp2;
//   else
//      p_lp = 0.5 * pre_left * (vel_left + fabs(vel_left)) / (vel_left);
//
////......pre_right-
//   temp1 = vel_right / c_mean - 1.0;
//   temp2 = 2.0 + vel_right / c_mean;
//   if(amach_r <= 1.0)
//      p_rm = 0.25 * pre_right * temp1 * temp1 * temp2;
//   else
//      p_rm = 0.5 * pre_right * (vel_right - fabs(vel_right)) / (vel_right);
//
//   ak_half = 10.0;
//   p_half = p_lp + p_rm;
//   p_min = (pre_left < pre_right) ? pre_left : pre_right;
//   s_half = fabs(pre_right - pre_left) / p_min;
//   a_min = ak_half * s_half;
//   a_min = (a_min < 1.0) ? a_min : 1.0;
//   s_half = 0.5 * a_min;
//
////......calculating the inviscid interface fluxes
//
////......mass flux
//   fl[0] = u_lp * rho_left + u_rm * rho_right;
//
////......enthalpy flux
//   fl[2] = 0.5 * (fl[0] * (h_left + h_right)
//                  - fabs(fl[0]) * (h_right - h_left));
//
////......momentum flux
//
//   flux_v1 = rho_left * vel_left * u_lp + rho_right * vel_right * u_rm;
//
//   flux_d1 = 0.5 * (fl[0] * (vel_left + vel_right) -
//                    fabs(fl[0]) * (vel_right - vel_left));
//
//   flux_p = p_half;
//
//   flux_n = (0.5 + s_half) * flux_v1 + (0.5 - s_half) * flux_d1;
//
//   fl[1] = flux_n + flux_p;
//}
