# Parallel Efficient Inner Product for Stabilizer States

**This project was originally submitted as a final project for the MasterMath course of Parallel Algorithms on January 11th, 2021. Corrections and improvements have been performed since then.** 

In this project we parallelize some of the algorithms related to the computation of the inner product for stabilizer states, a common computation when dealing with the classical simulation of quantum computers. The project is based on [this](https://arxiv.org/pdf/1210.6646.pdf) article by *H. J. García, I. G. Markov, and A. W. Cross*. There are a few differences between our implementation and theirs, as they consider an *nxn* matrix of Z literals, while we map each row of this matrix into a binary string of length *2n*. Every row is now an element of *GF(2)*, and so the arithmetic operations performed during the algorithm should be considered in terms of this field. This transformation is made with the intention of considering sparse implementations where the majority of the elements of this tableau are zeros, and the only non-zero elements are ones.
More information about the theory of stabilizer operators, methods and results can be found in the project report.

## Authors
* Ricardo Rivera Cardoso
* Andi Lin

## Software
The software used for this project is [**MulticoreBSP**](http://www.multicorebsp.com/download/) \[1\] and is complemented by the book *Parallel Scientific Computation: A Structured Approach Using BSP (2nd edition)* \[2\], which discusses the bases of parallel computation and applies these to many examples. Finally, [BSPedupack](https://webspace.science.uu.nl/~bisse101/Book2/psc2.html) provides beginner-friendly parallel algorithms that are implemented using this software.

\[1\] *A. N. Yzelman, R. H. Bisseling, D. Roose, K. Meerbergen, MulticoreBSP for C: a high-performance library for shared-memory parallel programming, technical report TW 624, revision May 2013, KU Leuven (post-print version).*
\[2\] *Bisseling, R. H. (2020). Parallel Scientific Computation: A Structured Approach Using BSP. Oxford University Press, USA.*

## License
This code and the MulticoreBSP software are both licensed under the GNU Lesser General Public License 3.0 - see [LICENSE.md](LICENSE.md) for more details.