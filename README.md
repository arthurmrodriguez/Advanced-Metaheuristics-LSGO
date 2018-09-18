# Advanced-Metaheuristics-LSGO

This repo contains my Computer Engineering Degree's Final Project, studied at the University of Granada, Spain. The main focus is to apply advanced metaheuristic algorithms into a **Big Optimization** problem with thousands of variables.
Our task is to find out how accurate are theoretical benchmark results compared to real EEG (Electroencephalography) data. 

This research is directed by Ph.D. Daniel Molina Cabrera.

The set of algorithms that will be compared are belong to a category known as Large Scale Global Optimization (LSGO), where the search space has more than a thousand of variables. These techniques have been designed to process theoretical benchmark functions. For more information, see the **Documentation** section to gain a broader view about this research.

Source code of each proposal has been obtained through either the director of this research or personal repositories of the authors.

### 1. MOS-SOCO2011
Multiple Offspring Sampling hybrid-based algorithm that combines a Differential Evolution Algorithm (DE) with a powerful Local Search (LS), the MTS-LS1. This algorithm was proposed for the SOCO 2011. See the following paper for further information: LaTorre, A., Muelas, S., Peña, J.-M. A MOS-based dynamic memetic differential evolution algorithm for continuous optimization: A scalability test (2011) Soft Computing, 15 (11), pp. 2187-2199

### 2. MOS-CEC2013
Enhanced version of the original MOS proposal, where a Classic Genetic Algorithm is combined along with a Solis Wets algorithm and a specifically designed Local Search for this proposal, the MTS-LS1-Reduced. This memetic algorithm has been considered since 2013 as the *state-of-the-art* algorithm within LSGO area. Paper of the proposal:
Latorre, A., Muelas, S., Pena, J.-M. Large scale global optimization: Experimental results with MOS-based hybrid algorithms
(2013) 2013 IEEE Congress on Evolutionary Computation, CEC 2013, art. no. 6557901, pp. 2742-2749

### 3. SHADEILS
SHADE-based memetic algorithm for LSGO where  one out of two powerful Local Search methods is selected according the needs of the optimization process: MTS-LS1 and L-BFGS-B. This algorithm is currently considered the *state-of-the-art* algorithms within LSGO. Paper of the publication: D. Molina, A. LaTorre, F. Herrera. SHADE with Iterative Local Search for Large-Scale Global Optimization. 2018 World Congress on Computational Intelligence (WCCI-2018), 2018 IEEE Conference on Evolutionary Computation (IEEE CEC'2018), Rio de Janeiro (Brasil), 1252-1259, July 8-13, 2018.

### 4. MLSHADE-SPA
Cooperative Co-evolution algorithm where three powerful Differential Evolution (DE) algorithms are used: LSHADE-SPA, EADE and ANDE. Furthermore, a modified version of the MTS algorithm is also used, MMTS. Paper of the publication (Under Publication) Anas A. Hadi, Ali W. Mohamed, and Kamal M. Jambi: LSHADE-SPA Memeteic Framework for Solving Large Scale Problems.

### 5. DG2
Differential Grouping algorithm for variable interaction detection. This proposal is an enhanced version of the original DG algorithm and aims to solve all the problems of the DG family algorithms. See the paper for further information: 
Omidvar, M.N., Yang, M., Mei, Y., Li, X., Yao, X. DG2: A faster and more accurate differential grouping for large-scale black-box optimization (2017) IEEE Transactions on Evolutionary Computation, 21 (6), pp. 929-942.













