
#Malandro

Malandro is a macro for the [CCPN Analysis software](http://www.ccpn.ac.uk/software/analysis) for the analysis of Nuclear Magnetic Resonance spectra. It is meant to help with sequential resonance assignment proteins. Sequential assignment of proteins in general consists of two steps:

1. grouping resonances that belong to the same residue together in units that are prefered to as _spin systems_ in CCPN Analysis, and 
2. determining the position of these spin system on the amino acid sequence of the molecule.

Especially the second step in the process can be very non-trivial and time consuming. This macro therefor focuses on helping with this step. I tried to make the algorithm as flexible as possible being able to use an arbitrary number of spectra with different magnetization transfer pathways, number of dimensions and labelling schemes. Therefor it can be useful in situations where there is a lot of information is present in a large number of different spectra that has to be combined in order to come up with a sequential assignment. The trivial idea is that the more peaks are supporting a sequential assignment and the better they fit the chemical shifts of the involved resonances, the more likely it is to be true. In a difficult assignment project it might not directly generate a full assignment. It might help to start of though, and once a few assignments have been made the algorithm can be run again in a repeated fashion. I wrote this macro with solid state NMR in the back of my mind, but don't see any reason why it would not come in handy in some cases in solution NMR. Of course, if there is an assignment algorithm that is tailor made for the type of spectra you have, I recommend you to use that.

##Introduction

I wrote this macro because I have a lot of different experiments and different labelling schemes. The basic thing I wanted to be able to quickly see is: If spin system A and spin system B are sequential neighbors, which peaks support this hypothesis in the different types of spectra and labelling schemes. If there is for instance a Threonine - Alanine pair in the sequence I wanted to find out which peaks I should observe in which spectra. Subsequently _all_ combinations of _each_ spin system that could be a Threonine and _each_ spin system that could be an Alanine spin should be made. The resonance frequencies in the spin systems can be used to check in the peak lists of the spectra to see which of these peaks indeed show up. The more peaks show up for a certain combination of spin systems, the more likely they are to be neighbors.
From here it was only a small step to extend this very local solution to the entire sequence and do a Monte Carlo / Annealing procedure to find overall good solutions.

I wrote the algorithm in such a way that it takes almost all necessary information from the analysis project. This prevents the in- and output a lot of data and makes it easy to run it after some new assignment are made. Seen from a high level the procedure follows the following steps (in there might be some steps in between to speed things up a bit):

1.  The expGraphs and Labelling Schemes of a set of selected experiments in used to calculate which correlations between atoms give rise to intra-residual and  sequential peaks in the sequence. The calculated co-labelling takes into account all labelled fractions of all isotopomers in every step of the magnetization transfer pathway.  

2. For each residue pair in the sequence, the spin systems that could fit the residue types are taken and all combinations are made to search for peaks in the peak lists that would fit within the tolerance of the peak dimension as configured in your _CCPN Analysis_ project. How well the peak fits is scored by a flat-bottomed potential, furthermore the symmetry of the spectrum is taken into account: for instance in a through space carbon-carbon correlation a peak containing the same information shows up on both sides of the diagonal. For spectra like these the score is divided by the number of peaks that have the same information, in a the carbon-carbon through-space correlation this is 2. Also the higher the dimensionality of the spectrum is, the higher the peak will be scored, this is done by multiplying with the square of the number of involved resonances. It would maybe be simpler to use the amount of dimensions. However, this does not reflect the fact that for instance a peak in a 3D spectrum like a NCOCX spectrum, a N-CO-CB cross peak is far more valuable than the partly diagonal N-CO-CO peak. Furthermore the overall probability that two spin systems are a sequencial pair is scored by the amount of resonances that contributes to all connecting peaks. The reasoning behind this last score is that two spin systems that potentially have a lot of connecting peaks that all have dimensional contributions from the same small set of resonances is not as good as a two spin systems that have the same amount, or even a slightly lower amount, of connecting peaks that have a larger set of different resonaces contributing to peak dimensions. 

3. The information that is generated in the first two steps is now used to perform a Monte Carlo / simulated annealing procedure. First a random assignment is made. This random assignment is completely consistent in the sense that there are no mismatched between the type of each residue in the sequence and the residue types a spin system could have. The assignment will stay consistent through-out the optimization. If there are less spin systems of a certain type than residues of this type _Joker_ spin systems will be inserted so that each residue will always have a spin system assigned to it. Spin systems however can have no residue assignment at all, in this case they are just waiting around until they are switched around with a spin system that _is_ assigned to a residue. In every step of the Monte Carlo procedure an attempt is made to switch two spin systems around. Before the optimization starts a set of 'allowed residues' is generated for each spin system based on its possible amino acid types and tentative assignments. Two spin systems can potentially be exchanged when there is an intersection between the 'allowed residues'sets. In the Monte Carlo procedure a random spin system is picked and a second one is randomly picked from a list of potential exchange partners. During the procedure it is checked whether this exchange is really valid considering the current assignment state of the two spin systems, to keep the solution consistent. The difference in score is calculated by summing over all relevant peak scores calculated in step 2 divided by the degeneracy of the peak. The degeneracy of a peak is simply defined as the amount of different contributions the peak has at the specific state of the assignment (depending on the peak, before or after the switch). Peaks can have theoretically have infinite degeneracy, but to prevent that these peaks influence the outcome of the optimization in a disproportionate way, their score is normalized in the described fashion. As in all Monte Carlo procedures the score is used to decide whether the switch is accepted or not: if the new score is better than the old, the switch is accepted and if it is worse it is accepted with a chance of e**(-delta score / kB*T). In practice I defined (at least for now) the value 1/(kB*T) as 'temperature constant', which can be manipulated in the temperature constant list in the GUI.



## Getting Started

### Installation

Compiled releases for linux and mac coming soon. For now, cython and a c-compiler, like gcc, have to be installed.

Cythonize and compile by running:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.bash}

python setup.py build_ext --inplace

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now the macro can be opened in CCPN analysis as usual by loading the start_malandro.py macro.


### The Random Generator

The included random generator is a re-implementation by Josh Ayers of the Mersenne Twister that is present in the normal python distribution. It is distributed as part of a library with other utilities for Cython: https://bitbucket.org/joshayers/cythonlib .

