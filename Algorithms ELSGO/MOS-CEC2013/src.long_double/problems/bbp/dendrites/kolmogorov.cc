#include <string>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>

#include <iostream>
#include <vector>

long double KolmogorovProb(long double z)
{
   // Calculates the Kolmogorov distribution function,
   //Begin_Html
   /*
   <img src="gif/kolmogorov.gif">
   */
   //End_Html
   // which gives the probability that Kolmogorov's test statistic will exceed
   // the value z assuming the null hypothesis. This gives a very powerful
   // test for comparing two one-dimensional distributions.
   // see, for example, Eadie et al, "statistocal Methods in Experimental
   // Physics', pp 269-270).
   //
   // This function returns the confidence level for the null hypothesis, where:
   //   z = dn*sqrt(n), and
   //   dn  is the maximum deviation between a hypothetical distribution
   //       function and an experimental distribution with
   //   n    events
   //
   // NOTE: To compare two experimental distributions with m and n events,
   //       use z = sqrt(m*n/(m+n))*dn
   //
   // Accuracy: The function is far too accurate for any imaginable application.
   //           Probabilities less than 10^-15 are returned as zero.
   //           However, remember that the formula is only valid for "large" n.
   // Theta function inversion formula is used for z <= 1
   //
   // This function was translated by Rene Brun from PROBKL in CERNLIB.

   long double fj[4] = {-2,-8,-18,-32}, r[4];
   const long double w = 2.50662827;
   // c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
   const long double c1 = -1.2337005501361697;
   const long double c2 = -11.103304951225528;
   const long double c3 = -30.842513753404244;

   long double u = fabs(z);//abs(z);
   long double p;
   if (u < 0.2) {
      p = 1;
   } else if (u < 0.755) {
      long double v = 1./(u*u);
      p = 1 - w*(exp(c1*v) + exp(c2*v) + exp(c3*v))/u;
   } else if (u < 6.8116) {
      r[1] = 0;
      r[2] = 0;
      r[3] = 0;
      long double v = u*u;
      int maxj = std::max(1, (int)round(3./u));
      for (int j=0; j<maxj;j++) {
         r[j] = exp(fj[j]*v);
      }
      p = 2*(r[0] - r[1] +r[2] - r[3]);
   } else {
      p = 0;
   }
   return p;
   }

//______________________________________________________________________________
long double KolmogorovTest(const std::vector<long double>& a, const std::vector<long double>& b, const char *option)
{
//  Statistical test whether two one-dimensional sets of points are compatible
//  with coming from the same parent distribution, using the Kolmogorov test.
//  That is, it is used to compare two experimental distributions of unbinned data.
//
//  Input:
//  a,b: One-dimensional arrays of length na, nb, respectively.
//       The elements of a and b must be given in ascending order.
//  option is a character string to specify options
//         "D" Put out a line of "Debug" printout
//         "M" Return the Maximum Kolmogorov distance instead of prob
//
//  Output:
// The returned value prob is a calculated confidence level which gives a
// statistical test for compatibility of a and b.
// Values of prob close to zero are taken as indicating a small probability
// of compatibility. For two point sets drawn randomly from the same parent
// distribution, the value of prob should be uniformly distributed between
// zero and one.
//   in case of error the function return -1
//   If the 2 sets have a different number of points, the minimum of
//   the two sets is used.
//
// Method:
// The Kolmogorov test is used. The test statistic is the maximum deviation
// between the two integrated distribution functions, multiplied by the
// normalizing factor (rdmax*sqrt(na*nb/(na+nb)).
//
//  Code adapted by Rene Brun from CERNLIB routine TKOLMO (Fred James)
//   (W.T. Eadie, D. Drijard, F.E. James, M. Roos and B. Sadoulet,
//      Statistical Methods in Experimental Physics, (North-Holland,
//      Amsterdam 1971) 269-271)
//
//  Method Improvement by Jason A Detwiler (JADetwiler@lbl.gov)
//  -----------------------------------------------------------
//   The nuts-and-bolts of the TMath::KolmogorovTest() algorithm is a for-loop
//   over the two sorted arrays a and b representing empirical distribution
//   functions. The for-loop handles 3 cases: when the next points to be
//   evaluated satisfy a>b, a<b, or a=b:
//
//      for (int i=0;i<na+nb;i++) {
//         if (a[ia-1] < b[ib-1]) {
//            rdiff -= sa;
//            ia++;
//            if (ia > na) {ok = true; break;}
//         } else if (a[ia-1] > b[ib-1]) {
//            rdiff += sb;
//            ib++;
//            if (ib > nb) {ok = true; break;}
//         } else {
//            rdiff += sb - sa;
//            ia++;
//            ib++;
//            if (ia > na) {ok = true; break;}
//            if (ib > nb) {ok = true; break;}
//        }
//         rdmax = max(rdmax,abs(rdiff));
//      }
//
//   For the last case, a=b, the algorithm advances each array by one index in an
//   attempt to move through the equality. However, this is incorrect when one or
//   the other of a or b (or both) have a repeated value, call it x. For the KS
//   statistic to be computed properly, rdiff needs to be calculated after all of
//   the a and b at x have been tallied (this is due to the definition of the
//   empirical distribution function; another way to convince yourself that the
//   old CERNLIB method is wrong is that it implies that the function defined as the
//   difference between a and b is multi-valued at x -- besides being ugly, this
//   would invalidate Kolmogorov's theorem).
//
//   The solution is to just add while-loops into the equality-case handling to
//   perform the tally:
//
//         } else {
//            long double x = a[ia-1];
//            while(a[ia-1] == x && ia <= na) {
//              rdiff -= sa;
//              ia++;
//            }
//            while(b[ib-1] == x && ib <= nb) {
//              rdiff += sb;
//              ib++;
//            }
//            if (ia > na) {ok = true; break;}
//            if (ib > nb) {ok = true; break;}
//         }
//
//  NOTE1
//  A good description of the Kolmogorov test can be seen at:
//    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm

   int na = a.size();
   int nb = b.size();

   std::string opt = option;
   std::transform(opt.begin(), opt.end(),opt.begin(), ::toupper);
   //opt.ToUpper();

   long double prob = -1;
//      Require at least two points in each graph
   if (na <= 2 || nb <= 2) {
      std::cerr << "KolmogorovTest: " << "Sets must have more than 2 points" << std::endl;
      return prob;
   }
//     Constants needed
   long double rna = na;
   long double rnb = nb;
   long double sa  = 1./rna;
   long double sb  = 1./rnb;
   long double rdiff;
   int ia,ib;
//     Starting values for main loop
   if (a[0] < b[0]) {
      rdiff = -sa;
      ia = 2;
      ib = 1;
   } else {
      rdiff = sb;
      ib = 2;
      ia = 1;
   }
   long double rdmax = fabs(rdiff);
//    Main loop over point sets to find max distance
//    rdiff is the running difference, and rdmax the max.
   bool ok = false;
   for (int i=0;i<na+nb;i++) {
      if (a[ia-1] < b[ib-1]) {
         rdiff -= sa;
         ia++;
         if (ia > na) {ok = true; break;}
      } else if (a[ia-1] > b[ib-1]) {
         rdiff += sb;
         ib++;
         if (ib > nb) {ok = true; break;}
      } else {
         long double x = a[ia-1];
         while(a[ia-1] == x && ia <= na) {
            rdiff -= sa;
            ia++;
         }
         while(b[ib-1] == x && ib <= nb) {
            rdiff += sb;
            ib++;
         }
         if (ia > na) {ok = true; break;}
         if (ib > nb) {ok = true; break;}
      }
      rdmax = std::max(rdmax,fabs(rdiff));
   }
//    Should never terminate this loop with ok = false!

   if (ok) {
      rdmax = std::max(rdmax,fabs(rdiff));
      long double z = rdmax * sqrt(rna*rnb/(rna+rnb));
      prob = KolmogorovProb(z);
   }
      // debug printout
   if ( opt.find("D")  != std::string::npos ) {
      printf(" Kolmogorov Probability = %g, Max Dist = %g\n",prob,rdmax);
   }
   if( opt.find("M")  != std::string::npos ) return rdmax;
   else                  return prob;
}
