/**
  * @file
  * @brief Common genetic operators.
  *
  */

#ifndef GACommonOps_H
#define GACommonOps_H

#include <vector>
#include <algorithm>

#include "genomes/GAMask.h"
#include "genomes/GA1DArrayGenome.h"

// Initializers

// The random initializer sets the elements of the array based on the alleles
// set.  We choose randomly the allele for each element.
template <typename T> void UniformInitializer (GAGenome& c) {

   GA1DArrayAlleleGenome<T>& child = DYN_CAST (GA1DArrayAlleleGenome<T>&, c);

   for (int i = child.length () - 1; i >= 0; i--)
      child.gene (i, child.alleleset (i).allele ());

}


// Random initializer for order-based genome.  Loop through the genome
// and assign each element the next allele in the allele set.  Once each
// element has been initialized, scramble the contents by swapping elements.
// This assumes that there is only one allele set for the array.
template <typename T> void OrderedInitializer (GAGenome& c) {

   GA1DArrayAlleleGenome<T>& child = DYN_CAST (GA1DArrayAlleleGenome<T>&, c);

   int length = child.length () - 1;
   int n = 0;
   int i;

   for (i = length; i >= 0; i--) {

      child.gene (i, child.alleleset ().allele (n++));

      if (n >= child.alleleset ().size ())
         n = 0;

   }

   for (i = length; i >= 0; i--)
      child.swap (i, GARandomInt (0, length));

}


// Comparators

// The comparator is supposed to return a number that indicates how similar
// two genomes are, so here we just compare elements and return a number that
// indicates how many elements match.  If they are different lengths then we
// return -1 to indicate that we could not calculate the differences.
// This assumes that there is an operator == defined for the object in the
// elements of the array.
template <typename T> long double ElementComparator (const GAGenome& a, const GAGenome& b) {

   const GA1DArrayGenome<T>& sis = DYN_CAST (const GA1DArrayGenome<T>&, a);
   const GA1DArrayGenome<T>& bro = DYN_CAST (const GA1DArrayGenome<T>&, b);

   if (sis.length () != bro.length ())
      return -1;

   if (sis.length () == 0)
      return 0;

   long double count = 0.0;

   for (int i = sis.length () - 1; i >= 0; i--)
      count += ((sis.gene (i) == bro.gene (i)) ? 0 : 1);

   return count / sis.length ();

}


// Crossovers

// Single point crossover for 1D array genomes.  Pick a single point then
// copy genetic material from each parent.  We must allow for resizable genomes
// so be sure to check the behaviours before we do the crossovers.  If resizing
// is allowed then the children will change depending on where the site is
// located.  It is also possible to have a mixture of resize behaviours, but
// we won't worry about that at this point.  If this happens we just say that
// we cannot handle that and post an error message.
template <typename T> int OnePointCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   unsigned int momsite, momlen;
   unsigned int dadsite, dadlen;

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE &&
          bro.resizeBehaviour () == GAGenome::FIXED_SIZE    ) {

         if (mom.length () != dad.length () ||
             sis.length () != bro.length () ||
             sis.length () != mom.length ()    ) {

            GAErr (GA_LOC, mom.className (), "one-point cross", gaErrSameLengthReqd);

            return nc;

         }

         momsite = dadsite = GARandomInt (0, mom.length ());
         momlen  = dadlen  = mom.length () - momsite;

      }
      else if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE ||
               bro.resizeBehaviour () == GAGenome::FIXED_SIZE    ) {

         GAErr (GA_LOC, mom.className (), "one-point cross", gaErrSameBehavReqd);

         return nc;

      }
      else {

         momsite = GARandomInt (0, mom.length ());
         dadsite = GARandomInt (0, dad.length ());

         momlen = mom.length () - momsite;
         dadlen = dad.length () - dadsite;

         sis.resize (momsite + dadlen);
         bro.resize (dadsite + momlen);

      }

      sis.copy (mom, 0, 0, momsite);
      sis.copy (dad, momsite, dadsite, dadlen);
      bro.copy (dad, 0, 0, dadsite);
      bro.copy (mom, dadsite, momsite, momlen);

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE) {

         if (mom.length () != dad.length () || sis.length () != mom.length ()) {

            GAErr (GA_LOC, mom.className (), "one-point cross", gaErrSameLengthReqd);
            return nc;

         }

         momsite = dadsite = GARandomInt (0, mom.length ());
         momlen  = dadlen  = mom.length () - momsite;

      }
      else {

         momsite = GARandomInt (0, mom.length ());
         dadsite = GARandomInt (0, dad.length ());

         momlen = mom.length () - momsite;
         dadlen = dad.length () - dadsite;

         sis.resize (momsite + dadlen);

      }

      if (GARandomBit ()) {

         sis.copy (mom, 0, 0, momsite);
         sis.copy (dad, momsite, dadsite, dadlen);

      }
      else {

         sis.copy (dad, 0, 0, dadsite);
         sis.copy (mom, dadsite, momsite, momlen);

      }

      nc = 1;

   }

   return nc;

}


// Two points crossover for the 1D array genome.  Similar to the single point
// crossover, but here we pick two points then grab the sections based upon
// those two points.
//   When we pick the points, it doesn't matter where they fall (one is not
// dependent upon the other).  Make sure we get the lesser one into the first
// position of our site array.
template <typename T> int TwoPointsCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   unsigned int momsite [2], momlen [2];
   unsigned int dadsite [2], dadlen [2];

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE &&
          bro.resizeBehaviour () == GAGenome::FIXED_SIZE    ) {

         if (mom.length () != dad.length () ||
             sis.length () != bro.length () ||
             sis.length () != mom.length ()    ) {

            GAErr (GA_LOC, mom.className (), "two-point cross", gaErrSameLengthReqd);
            return nc;

         }

         momsite [0] = GARandomInt (0, mom.length ());
         momsite [1] = GARandomInt (0, mom.length ());

         if (momsite [0] > momsite [1])
            std::swap (momsite [0], momsite [1]);

         momlen [0] = momsite [1] - momsite [0];
         momlen [1] = mom.length () - momsite [1];

         dadsite [0] = momsite [0];
         dadsite [1] = momsite [1];

         dadlen [0] = momlen [0];
         dadlen [1] = momlen [1];

      }
      else if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE ||
               bro.resizeBehaviour () == GAGenome::FIXED_SIZE    ) {

         return nc;

      }
      else {

         momsite [0] = GARandomInt (0, mom.length ());
         momsite [1] = GARandomInt (0, mom.length ());

         if (momsite [0] > momsite [1])
            std::swap (momsite [0], momsite [1]);

         momlen [0] = momsite [1] - momsite [0];
         momlen [1] = mom.length () - momsite [1];

         dadsite [0] = GARandomInt (0, dad.length ());
         dadsite [1] = GARandomInt (0, dad.length ());

         if (dadsite [0] > dadsite [1])
            std::swap (dadsite [0], dadsite [1]);

         dadlen [0] = dadsite [1] - dadsite [0];
         dadlen [1] = dad.length () - dadsite [1];

         sis.resize (momsite [0] + dadlen [0] + momlen [1]);
         bro.resize (dadsite [0] + momlen [0] + dadlen [1]);

      }

      sis.copy (mom, 0, 0, momsite [0]);
      sis.copy (dad, momsite [0], dadsite [0], dadlen [0]);
      sis.copy (mom, momsite [0] + dadlen [0], momsite [1], momlen [1]);

      bro.copy (dad, 0, 0, dadsite [0]);
      bro.copy (mom, dadsite [0], momsite [0], momlen [0]);
      bro.copy (dad, dadsite [0] + momlen [0], dadsite [1], dadlen [1]);

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      if (sis.resizeBehaviour () == GAGenome::FIXED_SIZE) {

         if (mom.length () != dad.length () || sis.length () != mom.length ()) {

            GAErr (GA_LOC, mom.className (), "two-point cross", gaErrSameLengthReqd);
            return nc;

         }

         momsite [0] = GARandomInt (0, mom.length ());
         momsite [1] = GARandomInt (0, mom.length ());

         if (momsite [0] > momsite [1])
            std::swap (momsite [0], momsite [1]);

         momlen [0] = momsite [1] - momsite [0];
         momlen [1] = mom.length () - momsite [1];

         dadsite [0] = momsite [0];
         dadsite [1] = momsite [1];

         dadlen [0] = momlen [0];
         dadlen [1] = momlen [1];

      }
      else {

         momsite [0] = GARandomInt (0, mom.length ());
         momsite [1] = GARandomInt (0, mom.length ());

         if (momsite [0] > momsite [1])
            std::swap (momsite [0], momsite [1]);

         momlen [0] = momsite [1] - momsite [0];
         momlen [1] = mom.length () - momsite [1];

         dadsite [0] = GARandomInt (0, dad.length ());
         dadsite [1] = GARandomInt (0, dad.length ());

         if (dadsite [0] > dadsite [1])
            std::swap (dadsite [0], dadsite [1]);

         dadlen [0] = dadsite [1] - dadsite [0];
         dadlen [1] = dad.length () - dadsite [1];

         sis.resize (momsite [0] + dadlen [0] + momlen [1]);

      }

      if (GARandomBit ()) {

         sis.copy (mom, 0, 0, momsite [0]);
         sis.copy (dad, momsite [0], dadsite [0], dadlen [0]);
         sis.copy (mom, momsite [0] + dadlen [0], momsite [1], momlen [1]);

      }
      else {

         sis.copy (dad, 0, 0, dadsite [0]);
         sis.copy (mom, dadsite [0], momsite [0], momlen [0]);
         sis.copy (dad, dadsite [0] + momlen [0], dadsite [1], dadlen [1]);

      }

      nc = 1;

   }

   return nc;

}


// Even and odd crossover for the array works just like it does for the
// binary strings.  For even crossover we take the 0th element and every other
// one after that from the mother.  The 1st and every other come from the
// father.  For odd crossover, we do just the opposite.
template <typename T> int EvenOddCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   int i;

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      if (sis.length () == bro.length () &&
          mom.length () == dad.length () &&
          sis.length () == mom.length ()    ) {

         for (i = sis.length () - 1; i >= 1; i -= 2) {

            sis.gene (i, mom.gene (i));
            bro.gene (i, dad.gene (i));
            sis.gene (i - 1, dad.gene (i - 1));
            bro.gene (i - 1, mom.gene (i - 1));

         }

         if (i == 0) {

            sis.gene (0, mom.gene (0));
            bro.gene (0, dad.gene (0));

         }

      }
      else {

         int start;
         int min = (mom.length () < dad.length ()) ? mom.length () : dad.length ();

         start = (sis.length () < min) ? sis.length () - 1 : min - 1;

         for (i = start; i >= 0; i--)
            sis.gene (i, ((i % 2 == 0) ? mom.gene (i) : dad.gene (i)));

         start = (bro.length () < min) ? bro.length () - 1 : min - 1;

         for (i = start; i >= 0; i--)
            bro.gene (i, ((i % 2 == 0) ? dad.gene (i) : mom.gene (i)));

      }

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      if (mom.length () == dad.length () && sis.length () == mom.length ()) {

         for (i = sis.length () - 1; i >= 1; i -= 2) {

            sis.gene (i, mom.gene (i));
            sis.gene (i - 1, dad.gene (i - 1));

         }

         if (i == 0)
            sis.gene (0, mom.gene (0));

      }
      else {

         int min = (mom.length () < dad.length ()) ? mom.length () : dad.length ();
         min = (sis.length () < min) ? sis.length () - 1 : min - 1;

         for (i = min; i >= 0; i--)
            sis.gene (i, ((i % 2 == 0) ? mom.gene (i) : dad.gene (i)));

      }

      nc = 1;

   }

   return nc;

}


// Partial match crossover for the 1D array genome.  This uses the partial
// matching algorithm described in Goldberg's book.
//   Parents and children must be the same size for this crossover to work.  If
// they are not, we post an error message.
//   We make sure that b will be greater than a.
template <typename T> int PartialMatchCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   int a = GARandomInt (0, mom.length ());
   int b = GARandomInt (0, dad.length ());

   if (b < a)
      std::swap (a, b);

   int i, j, index;

   if (mom.length () != dad.length ()) {

      GAErr (GA_LOC, mom.className (), "parial match cross", gaErrBadParentLength);
      return nc;

   }

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      sis.GAArray<T>::copy (mom);

      for (i = a, index = a; i < b; i++, index++) {

         for (j = 0; j < sis.length () - 1 && sis.gene (j) != dad.gene (index); j++);
         sis.swap (i, j);

      }

      bro.GAArray<T>::copy (dad);

      for (i = a, index = a; i < b; i++, index++) {

         for (j = 0; j < bro.length () - 1 && bro.gene (j) != mom.gene (index); j++);
         bro.swap (i, j);

      }

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      const GA1DArrayGenome<T> *parent1, *parent2;

      if (GARandomBit ()) {

         parent1 = &mom;
         parent2 = &dad;

      }
      else {

         parent1 = &dad;
         parent2 = &mom;

      }

      sis.GAArray<T>::copy (*parent1);

      for (i = a, index = a; i < b; i++, index++) {

         for (j = 0; j < sis.length () - 1 && sis.gene (j) != parent2->gene (index); j++);
         sis.swap (i, j);

      }

      nc = 1;

   }

   return nc;

}


// Randomly take bits from each parent.  For each bit we flip a coin to see if
// that bit should come from the mother or the father.  If strings are
// different lengths then we need to use the mask to get things right.
template <typename T> int UniformCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int n = 0;
   int i;

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      if (sis.length () == bro.length () &&
          mom.length () == dad.length () &&
          sis.length () == mom.length ()    ) {

            for(i = sis.length () - 1; i >= 0; i--)
               if (GARandomBit ()) {

                  sis.gene (i, mom.gene (i));
                  bro.gene (i, dad.gene (i));

               }
               else {

                  sis.gene (i, dad.gene (i));
                  bro.gene (i, mom.gene (i));

               }

      }
      else {

         GAMask mask;

         int start;
         int max = (sis.length () > bro.length ()) ? sis.length () : bro.length ();
         int min = (mom.length () < dad.length ()) ? mom.length () : dad.length ();

         mask.size (max);

         for (i = 0; i < max; i++)
            mask [i] = GARandomBit ();

         start = (sis.length () < min) ? sis.length () - 1 : min - 1;

         for (i = start; i >= 0; i--)
            sis.gene (i, (mask [i] ? mom.gene (i) : dad.gene (i)));

         start = (bro.length () < min) ? bro.length () - 1 : min - 1;

         for (i = start; i >= 0; i--)
            bro.gene (i, (mask [i] ? dad.gene (i) : mom.gene (i)));

      }

      n = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      if (mom.length () == dad.length () && sis.length () == mom.length ())
         for (i = sis.length () - 1; i >= 0; i--)
            sis.gene (i, (GARandomBit () ? mom.gene (i) : dad.gene (i)));

      else {

         int min = (mom.length () < dad.length ()) ? mom.length () : dad.length ();

         min = (sis.length () < min) ? sis.length () : min;

         for (i = min - 1; i >= 0; i--)
            sis.gene (i, (GARandomBit () ? mom.gene (i) : dad.gene (i)));

      }

      n = 1;

   }

  return n;

}


// This function determines whether or not an indexed position is a hole that
// we can substitute into.  It does a linear search to find the holes (yuk).
template <typename T> int GA1DArrayIsHole (const GA1DArrayGenome<T>& c, const GA1DArrayGenome<T> &dad,
                                           int index, int a, int b) {

   for (int i = a; i < b; i++)
      if (c.gene (index) == dad.gene (i))
         return 1;

   return 0;

}


// Order crossover for the 1D array genome.  This uses the order crossover
// described in Goldberg's book.
//   Parents and children must be the same length.
//   We make sure that b will be greater than a.
//   This implementation isn't terribly smart.  For example, I do a linear
// search rather than caching and doing binary search or smarter hash tables.
//   First we copy the mother into the sister.  Then move the 'holes' into the
// crossover section and maintain the ordering of the non-hole elements.
// Finally, put the 'holes' in the proper order within the crossover section.
// After we have done the sister, we do the brother.
template <typename T> int OrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   int a = GARandomInt (0, mom.length ());
   int b = GARandomInt (0, mom.length ());

   if (b < a)
      std::swap (a, b);

   int i, j, index;

   if (mom.length () != dad.length ()) {

      GAErr (GA_LOC, mom.className (), "order cross", gaErrBadParentLength);
      return nc;

   }

   if (c1 && c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      // Copy the parent
      sis.GAArray<T>::copy (mom);

      // Move all the 'holes' into the crossover section
      for (i = 0, index = b; i < sis.size (); i++, index++) {

         if (index >= sis.size ())
            index = 0;

         if (GA1DArrayIsHole (sis, dad, index, a, b))
            break;

      }

      for (; i < sis.size () - b + a; i++, index++) {

         if (index >= sis.size ())
            index = 0;

         j = index;

         do {

            j++;

            if (j >= sis.size ())
               j = 0;

         } while (GA1DArrayIsHole (sis, dad, j, a, b));

         sis.swap (index, j);

      }

      // Now put the 'holes' in the proper order within the crossover section.
      for (i = a; i < b; i++)
         if (sis.gene (i) != dad.gene (i))
            for (j = i + 1; j < b; j++)
               if (sis.gene (j) == dad.gene (i))
                  sis.swap (i, j);

      // Now do the other child
      bro.GAArray<T>::copy (dad);

      // Move all the 'holes' into the crossover section
      for (i = 0, index = b; i < bro.size (); i++, index++) {

         if (index >= bro.size ())
            index = 0;

         if (GA1DArrayIsHole (bro, mom, index, a, b))
            break;

      }

      for(; i<bro.size()-b+a; i++, index++) {

         if (index >= bro.size ())
            index = 0;

         j = index;

         do {

            j++;

            if (j >= bro.size ())
               j = 0;

         } while (GA1DArrayIsHole (bro, mom, j, a, b));

         bro.swap (index, j);

      }

      // Now put the 'holes' in the proper order within the crossover section.
      for (i = a; i < b; i++)
         if (bro.gene (i) != mom.gene (i))
            for (j = i + 1; j < b; j++)
               if (bro.gene (j) == mom.gene (i))
                  bro.swap (i, j);

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      const GA1DArrayGenome<T> *parent1, *parent2;

      if (GARandomBit ()) {

         parent1 = &mom;
         parent2 = &dad;

      }
      else {

         parent1 = &dad;
         parent2 = &mom;

      }

      sis.GAArray<T>::copy (*parent1);

      for (i = 0, index = b; i < sis.size (); i++, index++) {

         if (index >= sis.size ())
            index = 0;

         if (GA1DArrayIsHole (sis, *parent2, index, a, b))
            break;

      }

      for (; i < sis.size () - b + a; i++, index++) {

         if (index >= sis.size ())
            index = 0;

         j = index;

         do {

            j++;

            if (j >= sis.size ())
               j = 0;

         } while (GA1DArrayIsHole (sis, *parent2, j, a, b));

         sis.swap (index, j);
      }

      for (i = a; i < b; i++)
         if (sis.gene (i) != parent2->gene (i))
            for (j = i + 1; j < b; j++)
               if (sis.gene (j) == parent2->gene (i))
                  sis.swap (i, j);

      nc = 1;

   }

   return nc;

}


// Cycle crossover for the 1D array genome.  This is implemented as described
// in goldberg's book.  The first is picked from mom, then cycle using dad.
// Finally, fill in the gaps with the elements from dad.
//   We allocate space for a temporary array in this routine.  It never frees
// the memory that it uses, so you might want to re-think this if you're really
// memory-constrained (similar to what we do with the uniform crossover when
// the children are resizeable).
//  Allocate space for an array of flags.  We use this to keep track of whether
// the child's contents came from the mother or the father.  We don't free the
// space here, but it is not a memory leak.
//   The first step is to cycle through mom & dad to get the cyclic part of
// the crossover.  Then fill in the rest of the sis with dad's contents that
// we didn't use in the cycle.  Finally, do the same thing for the other child.
//   Notice that this implementation makes serious use of the operator= for the
// objects in the array.  It also requires the operator != and == comparators.
template <typename T> int CycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2){

   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int nc = 0;
   int i, current = 0;

   if (mom.length () != dad.length ()) {

      GAErr (GA_LOC, mom.className (), "cycle cross", gaErrBadParentLength);
      return nc;

  }

   if (c1 && c2) {

      GAMask mask;
      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c1);
      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      mask.size (sis.length ());
      mask.clear ();

      sis.gene (0, mom.gene (0));
      mask [0] = 1;

      while (dad.gene (current) != mom.gene (0))
         for (i = 0; i < sis.size (); i++)
            if (mom.gene (i) == dad.gene (current)) {

               sis.gene (i, mom.gene (i));
               mask [i] = 1;
               current = i;
               break;

            }

      for (i = 0; i < sis.size (); i++)
         if (mask [i] == 0)
            sis.gene (i, dad.gene (i));

      mask.clear ();

      bro.gene (0, dad.gene (0));
      mask [0] = 1;

      while (mom.gene (current) != dad.gene (0))
         for (i = 0; i < bro.size (); i++)
            if (dad.gene (i) == mom.gene (current)) {

               bro.gene (i, dad.gene (i));
               mask [i] = 1;
               current = i;
               break;

            }


      for (i = 0; i < bro.size (); i++)
         if (mask [i] == 0)
            bro.gene (i, mom.gene (i));

      nc = 2;

   }
   else if (c1 || c2) {

      GA1DArrayGenome<T>& sis = (c1 ? DYN_CAST (GA1DArrayGenome<T>&, *c1) :
                                      DYN_CAST (GA1DArrayGenome<T>&, *c2)   );

      const GA1DArrayGenome<T> *parent1, *parent2;

      if (GARandomBit ()) {

         parent1 = &mom;
         parent2 = &dad;

      }
      else {

         parent1 = &dad;
         parent2 = &mom;

      }

      GAMask mask;
      mask.size (sis.length ());
      mask.clear ();

      sis.gene (0, parent1->gene (0));
      mask [0] = 1;

      while (parent2->gene (current) != parent1->gene (0))
         for (i = 0; i < sis.size (); i++)
            if (parent1->gene (i) == parent2->gene (current)) {

               sis.gene (i, parent1->gene (i));
               mask [i] = 1;
               current = i;
               break;

            }

      for (i = 0; i < sis.size (); i++)
         if (mask [i] == 0)
            sis.gene (i, parent2->gene (i));

      nc = 1;

   }

   return nc;

}


template <typename T> int AlternativeCycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   int n = 0;

   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int len = dad.length ();
   int pos = GARandomInt (0, len - 1);

   if (c1) {

      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c1);

      for (int i = 0; i < len; i++)
         bro.gene (i, dad.gene (i));

      int j = pos;

      do {

         bro.gene (j, mom.gene (j));
         int i = 0;

         while (dad.gene (j) != mom.gene (i))
            i++;

         j = i;

      } while (j != pos);

      n++;

   }

   if (c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      for (int i = 0; i < len; i++)
         sis.gene (i, mom.gene (i));

      int j = pos;

      do {

         sis.gene (j, dad.gene (j));
         int i = 0;

         while (mom.gene (j) != dad.gene (i))
            i++;

         j = i;

      } while (j != pos);

      n++;

   }

   return n;

}


template <typename T> int AlternativeOrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

   int n = 0;

   const GA1DArrayGenome<T>& dad = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& mom = DYN_CAST (const GA1DArrayGenome<T>&, p2);

   int len  = dad.length ();
   int pos1 = GARandomInt (0, len - 1);
   int pos2 = pos1;

   while (pos1 == pos2)
      pos2 = GARandomInt (0, len - 1);

   if (pos2 < pos1)
      std::swap (pos1, pos2);

   if (c1) {

      GA1DArrayGenome<T>& bro = DYN_CAST (GA1DArrayGenome<T>&, *c1);

      for (int i = 0; i < len; i++)
         bro.gene (i, mom.gene (i));

      std::vector<T> seq;

      for (int i = pos1 + 1; i <= pos2; i++)
         seq.push_back (mom.gene (i));

      int j = (pos2 + 1) % len;
      int k = (pos2 + 1) % len;

      while (j != pos1 + 1) {

         while (find (seq.begin (), seq.end (), dad.gene (k)) != seq.end ())
            k = (k + 1) % len;

         bro.gene (j, dad.gene (k));

         k = (k + 1) % len;
         j = (j + 1) % len;

      }

      seq.clear ();
      n++;

   }

   if (c2) {

      GA1DArrayGenome<T>& sis = DYN_CAST (GA1DArrayGenome<T>&, *c2);

      for (int i = 0; i < len; i++)
         sis.gene (i, dad.gene (i));

      std::vector<T> seq;

      for (int i = pos1 + 1; i <= pos2; i++)
         seq.push_back (dad.gene (i));

      int j = (pos2 + 1) % len;
      int k = (pos2 + 1) % len;

      while (j != pos1 + 1) {

         while (find (seq.begin (), seq.end (), mom.gene (k)) != seq.end ())
            k = (k + 1) % len;

         sis.gene (j, mom.gene (k));

         k = (k + 1) % len;
         j = (j + 1) % len;

      }

      seq.clear ();
      n++;

   }

   return n;

}


// DE Crossovers

template <typename T> int ExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const long double prob) {

   const GA1DArrayGenome<T>& x_i    = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   GA1DArrayGenome<T>&       x_new  = DYN_CAST (GA1DArrayGenome<T>&, *c1);

   int dim   = x_i.length();
   int nrand = GARandomInt(0,dim-1);
   int lrand = 0;

   while (GARandomDouble(0,1) < prob and (lrand<dim) )
    lrand = lrand+1;
   if (lrand==0)
    lrand=1;

   int initial_pos = ( nrand + lrand ) % dim;
   while (initial_pos != nrand) {
      x_new.gene(initial_pos,x_i.gene(initial_pos) );
      initial_pos = (initial_pos + 1) % dim;
   }

   return 1;

}

template <typename T> int BinomialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const long double prob) {

   const GA1DArrayGenome<T>& x_i    = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   const GA1DArrayGenome<T>& x_newc = DYN_CAST (const GA1DArrayGenome<T>&, p2);
   GA1DArrayGenome<T>&       x_new  = DYN_CAST (GA1DArrayGenome<T>&, *c1);

   int dim   = x_i.length();
   int irand = GARandomInt(0,dim);

   for(int j=0; j<dim; j++){

      long double tmprand = GARandomDouble(0.0,1.0);

      if( tmprand <= prob or j==irand )
         x_new.gene( j , x_newc.gene(j) );
      else
         x_new.gene( j , x_i.gene(j)    );

   }

   return 1;

}

// Mutators

// Randomly swap elements in the array.
template <typename T> int SwapMutator (GAGenome& c, long double pmut) {

   GA1DArrayGenome<T>& child = DYN_CAST (GA1DArrayGenome<T>&, c);

   //register int n, i;

   if (pmut <= 0.0)
      return (0);

   long double nMut = pmut * STA_CAST (long double, child.length ());
   int length = child.length () - 1;

//   if (nMut < 1.0) { // we have to do a flip test on each bit
//
//      nMut = 0;
//
//      for (i = length; i >= 0; i--)
//         if (GAFlipCoin (pmut)) {
//
//            child.swap (i, GARandomInt (0, length));
//            nMut++;
//
//      }
//
//   }
//   else  // only flip the number of bits we need to flip
   if (GAFlipCoin(nMut - (int)nMut) ) nMut++;

  for (int n = 0; n < (int) nMut; n++) {
    child.swap (GARandomInt (0, length), GARandomInt (0, length));
  }

  return (STA_CAST (int, nMut));

}


// Randomly pick elements in the array then set the element to any of the
// alleles in the allele set for this genome.  This will work for any number
// of allele sets for a given array.
template <typename T> int FlipMutator (GAGenome& c, long double pmut) {

   GA1DArrayAlleleGenome<T>& child = DYN_CAST (GA1DArrayAlleleGenome<T>&, c);

   register int n, i;

   if (pmut <= 0.0)
      return(0);

   long double nMut = pmut * STA_CAST (long double, child.length ());

//    if (nMut < 1.0) { // we have to do a flip test on each bit

//       nMut = 0;

//       for (i = child.length () - 1; i >= 0; i--)
//          if (GAFlipCoin (pmut)) {

//             child.gene(i, child.alleleset(i).allele());
//             nMut++;

//          }

//    }
//    else // only flip the number of bits we need to flip
   if (GAFlipCoin(nMut - (int) nMut)) nMut++;

   for (n = 0; n < nMut; n++) {

         i = GARandomInt (0, child.length () - 1);
         child.gene (i, child.alleleset (i).allele ());

   }

   return (STA_CAST (int, nMut));

}


template <typename T> int RepeatedExchangeMutator (GAGenome& g, long double pmut) {

   int nMut = 0;

   if (pmut <= 0.0)
      return 0;

   GA1DArrayGenome<T>& gen = DYN_CAST (GA1DArrayGenome<T>&, g);

   int len = gen.length ();

   for (int i = 0; i < len; i++) {

      if (GAFlipCoin (pmut)) {

         gen.swap (i, GARandomInt (0, len - 1));

         nMut++;

      }

   }

   return nMut;

}


template <typename T> int SimpleInversionMutator (GAGenome& g, long double pmut) {

   int nMut = 0;

   if (pmut <= 0.0)
      return 0;

   GA1DArrayGenome<T>& gen = DYN_CAST (GA1DArrayGenome<T>&, g);

   if (GAFlipCoin (pmut)) {

      int len   = gen.length ();
      int begin = GARandomInt (0, len - 1);
      int end   = begin;

      while (end == begin)
         end = GARandomInt (0, len - 1);

      // TODO: convert this to an int array to increase performance
      std::vector<T> subsequence;

      int i = begin;

      do {

         subsequence.push_back (gen.gene (i));
         i = (i + 1) % len;

      } while (i != ((end + 1) % len));

      i = begin;

      do {

         gen.gene (i, subsequence.back ());
         subsequence.pop_back ();
         i = (i + 1) % len;

      } while (i != ((end + 1) % len));

      nMut++;

   }

   return nMut;

}


#endif
