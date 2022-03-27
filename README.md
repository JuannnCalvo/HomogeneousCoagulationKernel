All the codes in this repository have been developed under the framework of the COCOMAS research project. This research has been funded by the spanish MINECO (Ministerio de Ciencia e Innovación), grant ref. MTM2017- 91054-EXP. These codes are provided under the GNU General Public License v3.0.


These codes costitute a proof of concept for a collocation+Tikhonov regularization approach to the inverse problem for homogeneous kernels in the binary coagulation equation and some oceanographical models based on it; as such, the codes are not optimized. For systhematic usage a separate pre-computation of the regularization matrices is recommended.

The collocation+Tikhonov regularization strategy was introduced in [R. Muralidar, D. Ramkrishna, J. Coll.  Interface Sci. 112 (1986), 348] for the inversion of homogeneous coagulation kernels.

These codes are based on several specific functionalities of the GNU Scientific Library (https://www.gnu.org/software/gsl/); release GSL 2.6 was used to design and test the routines.

General parameters of the colocation method can be tuned in the header files. Most of the routines use as input the homogeneity degree of the kernel. This can be guessed in advance by fitting temporal data as explained for instance in [R. Muralidar, D. Ramkrishna, J. Coll.  Interface Sci. 112 (1986), 348].

The optimal regularization parameter lambda is chosen using both the L-curve and the GCV-curve criteria; both reconstructions are provided. The GCV-curve criterion seems to perform poorly in this inversion problem.
