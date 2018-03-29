/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "solisrand.h"
#include <string>

void SRandom::setSeed(unsigned long seed) {
   m_seed = seed;
   changed = 1;
}

#define M  714025
#define IA 1366
#define IC 150889

SRandom::SRandom(unsigned long seed) {
    setSeed(seed);
}

long double  SRandom::rand(void) {
 	int j;
 	static long ir[98],iy;// idum = 13171191;
	static long idum = m_seed;
 	static int iff=0;
 
 	if (idum < 0 || iff == 0 || changed == 1) {
 		iff=1;
 		if ((idum=(IC-(idum)) % M) < 0) idum = -(idum);
 		for (j=1;j<=97;j++) {
 			idum=(IA*(idum)+IC) % M;
 			ir[j]=(idum);
 		}
 		idum=(IA*(idum)+IC) % M;
 		iy=(idum);
		changed = 0;
 	}
 	j=(int)(1 + 97.0*(iy)/M);  
 
 	if (j > 97 || j < 1)
 		throw new std::string("Failure in random number generator");
 
 	iy=ir[j];
 	idum=(IA*(idum)+IC) % M;
 	ir[j]=(idum);
	m_seed = idum;
 
 	return (float) (iy)/M;
}

unsigned long  SRandom::getSeed(void) {
   return m_seed;
}
