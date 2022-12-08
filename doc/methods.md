# Numerical methods

This document briefly describes the features of the spatial discretizations
available in MPI-AMRVAC. The different options can be set in the
[methodlist](@ref par_methodlist) of the input par file. For a more
extensive description, you can read the article [Comparison of some FCT and
TVD Schemes](http://www-personal.umich.edu/~gtoth/Papers/vac.html). Also, the
paper using MPI-AMRVAC has info on the various methods, see 'Parallel, grid-
adaptive approaches for relativistic hydro and magnetohydrodynamics', R.
Keppens, Z. Meliani, A.J. van Marle, P. Delmont, A. Vlasis, &amp; B. van der
Holst, 2011, JCP. [doi:10.1016/j.jcp.2011.01.020](http://dx.doi.org/10.1016/j.jcp.2011.01.020).
A more recent update is 'MPI-AMRVAC FOR SOLAR AND ASTROPHYSICS', O. Porth, C.
Xia, T. Hendrix, S. P. Moschou and R. Keppens ApJS.
[doi:10.1088/0067-0049/214/1/4](http://dx.doi.org/10.1088/0067-0049/214/1/4).

The acronyms TVD, TVDLF, and TVD-MUSCL stand for Total Variation Diminishing,
TVD Lax-Friedrich, and TVD Monotonic Upwind Scheme for Conservation Laws,
respectively. Then, depending on the physics module selected, you also have
HLL and HLLC schemes, which are due to harten Lax, van Leer, with the HLLC
variant including a treatment for the Contact discontinuity, as e.g. described
for Euler gas dynamics in _E.F. Toro, Riemann solvers and numerical methods
for fluid dynamics (Berlin, Springer-Verlag, 1997)_.

Not all methods are available or meaningfull for all physics modules. In fact,
we have the following combinations typically:

    Physics   Schemes
    --------------------------------------------------------------------------
    rho       TVDLF, HLL, HLLC, TVD (Roe solver), TVDMU (Roe solver), FD
    hd        TVDLF, HLL, HLLC, TVD (Roe solver), TVDMU (Roe solver), FD
    mhd       TVDLF, HLL, HLLC, TVD (Roe solver), TVDMU (Roe solver), FD, HLLD

Also, the method can be selected per AMR grid level, but one can not combine
different stepsize methods (hence, TVD is the only second order onestep
method, while all others can be used with twostep or fourstep typeadvance
setting). In MPI-AMRVAC, the **flux_scheme** is thus an array of strings, one
string per level up to **nlevelshi**. Some more info follows on the various
methods.

## 2nd Order Central Difference Scheme: flux_scheme='cd',...

The explicit central differencing schemes are not stable by themselves for
advection dominated problems. The second order central difference scheme
('cd') is used by the TVD scheme before the limiting is applied. Otherwise it
is useful for testing a few time steps, since this scheme contains no
artificial fluxes, thus comparison with analytic formulae is straightforward.
It is straightforward to generalize this central difference approach to higher
order accuracy, at the expense of introducing a wider stencil.

## High order finite difference Scheme: flux_scheme='fd',...

This implements conservative finite differences with global Lax-Friedrich flux
splitting. It can be used with almost all limiters (exluding ppm) and yields
high order accuracy in space. For second, third and fifth order reconstruction
you can set e.g.: **limiter='koren'/'cada3'/'mp5'**.

## TVDLF Scheme: flux_scheme='tvdlf'...

The TVD Lax-Friedrich method is robust, in most cases there are no spurious
oscillations, but it is somewhat more diffusive than other TVD or HLLC
methods. Since it does not use a Riemann solver, it is usually faster than TVD
or TVD-MUSCL.

The Courant number should be less than 1, **courantpar=0.8** is recommended.
Second order time discretization is best achieved by a Hancock predictor step,
so the corresponding **typepred1='hancock'**.

TVDLF can be used with **dimsplit=F**, it is also preferred for steady state
calculations.

The second order TVDLF scheme **flux_scheme='tvdlf'** uses limiters. There are
many choices available: the 'minmod' limiter gives the smoothest result, the
'woodward' limiter is sharper, and the 'superbee' limiter is probably too
sharp. The **'woodward'** limiter is recommended, but note that the default is
the most robust **limiter='minmod'**. The various options can be found
in the `mod_finite_volume.t` module, in the subroutine `dwlimiter2`. The 
slope limiting is performed on the primitive variables. 
You can even employ limiting on logarithmically stretched
variables (which should be positive, like a density or pressure), by setting
the `loglimit` flags. You can also use third order accurate
**limiter='ppm'**, but the code will run with a wider ghost
cell region, namely **nghostcells=4**. A third order limiter without a need 
to widen the ghost cell layers is the _'cada3'_ limiter (sometimes called LIMO3).

## TVD-MUSCL Scheme: flux_scheme='tvdmu'...

The TVD-MUSCL scheme is a two-step TVD algorithm using the same Hancock
predictor step and upwinding as TVDLF, and a characteristic based Riemann
solver similar to the TVD method. At the moment Riemann solvers are
implemented for adiabatic hydrodynamics, hydrodynamics, and full MHD. 
The scalar transport equation has a trivial Riemann solver. The scheme has 
comparable resolution to the non-MUSCL TVD method.

The Courant number should be less than 1, **courantpar=0.8** is recommended.

TVD-MUSCL can be dimensionally split **dimsplit=T** or unsplit **dimsplit=F**. 
The multistep Runge-Kutta schemes can be applied, such as **time_integrator='fourstep'**.

Linear Riemann solvers can produce non-physical solutions. This can be
eliminated by the use of an entropy fix, controlled by **typeentropy** and the
**entropycoef**. The default is **typeentropy='nul'**. See the details for
the entropy fixes in the respective `mod_PHYS_roe.t` files, as well as
in the `mod_tvd.t` module.

## TVD Scheme: flux_scheme='tvd',...

The non-MUSCL TVD method with Roe approximate Riemann solver is one of the
most accurate and efficient of the implemented schemes.

There are a few variants of the TVD scheme, but the default is
**typetvd='roe'**. Details are in the _mod_tvd.t _ module.

This solver has to be dimensionally split, set **dimsplit=T**.

The Courant number should be less than 1, **courantpar=0.8** is recommended.

The same limiters can be used as for TVDLF and TVD MUSCL, but they are applied
to the characteristic waves, rather than to the primitive variables. The
order of the characteristic waves is defined in the **mod_PHYS_roe.t*** files.
The **'woodward'** limiter is recommended, but note that the default is
**limiter='minmod'**.

The entropy fix for the Riemann solver is given by the **typeentropy** array,
it has the same meaning as for the TVD-MUSCL method, and for MHD, the
divergence B problem should also be taken care of.

## HLL and HLLC schemes

The TVDLF scheme hence uses minimal info on the wave speeds, and in
combination with AMR and its inherent robustness due to its diffusive nature,
it is readily usable for any system of conservation laws at minimal
implementation costs. Maximal wave speed info is used in a full Roe-type
approximate Riemann solver as employed by TVD or TVD-MUSCL, where all
characteristic wave speeds (7 in total for (relativistic) MHD) as well as the
wave strengths are deduced from the eigenvalues, as well as right and left
eigenvector pairs of the flux Jacobian. The simpler HLL, HLLC solvers, make
further approximations to their corresponding representation of the Riemann
fan, as schematically illustrated below. ![](figmovdir/solvers.gif) These type
of solvers originated in gas dynamics and Newtonian MHD, and have meanwhile
been adapted to relativistic (M)HD. Depending on the amount of waves used to
approximate the actual 7-wave fan, a corresponding amount of different fluxes
are computed. One switches between their expressions according to the relative
orientation of the wave signals in the (x,t) diagram. Appropriate recipes for
computing meaningful intermediate states ensure desirable properties like
positivity (positive pressures and densities), the ability to capture isolated
discontinuities, etc. For most physics modules, these HLL and HLLC variants are
available too. The HLLD variant is only applicable for MHD.

## Maintaining the div B=0 Condition

In multidimensional MHD the numerical conservation of divergence of magnetic field
div B is not guaranteed by the standard TVD or HLL type schemes. This can lead to 
inaccuracies as well as instablilities. For all the schemes below, you can 
influence how to compute div B, by setting _typegrad_ and _typediv_, along with 
_gradient_limiter_.
This allows to select either a standard central difference evaluation, or one
evaluated after the cell-center values have been reconstructed to the cell
edges. User can select one of the following methods by select **typedivbfix**
and related parameters in _mhd_list_ of par file.

#### Powell fix: typedivbfix='powel'

For multidimensional MHD calculations the non-conservative form of the [MHD
equations](@ref eq_mhd) seems to produce better results than the usual
conservative form. The idea is to include source terms proportional to div B
into the momentum, energy and induction equations and to add a divergence wave
for the Riemann solver.

Powell scheme is fast, it stabilizes the Riemann solver, and improves
results for TVDLF and similar type methods, but it is non-conservative, and
div B is not kept close to zero. 

#### Janhunen fix: typedivbfix='janhunen'

Source term in Powell fix is only added to the induction equation. This approach
is usable for both classical and relativistic MHD.

#### Diffusive fix: typedivbfix='linde'

You can also use the diffusive (parabolic) approach, see the
[equations](@ref eq_divb_fix). It uses a `C_d` coefficient quantified by
`divbdiff`, which can be up to 2. This method is by default inactive,
identified by `divbdiff=1`, but it is recommended for many multi-D MHD
applications.

#### Dedner fix: typedivbfix='glm1', 'glm2', or 'glm3'

This implements the mixed hyperbolic and parabolic dampening of the divB error
using an additional scalar variable _Psi_ (need an addition of the name and
boundary condition type in your par-file). The algorithm is described by
Dedner et al. in _Journal of Computational Physics 175, 645-673 (2002)
doi:10.1006/jcph.2001.6961_. The three versions differ in the source terms 
taken along. Thus 'glm1' corresponds 
to _Equation (24)_ of Dedner et al and 'glm2'
corresponds to _Equation (38)_ of this paper. The option 'glm3' adds no
additional sources to the MHD system. We recommend the option
'glm1'. For example: in your par-file,

    &mhd_list
    typedivbfix='glm1'
    ...

in your `mod_usr.t`, add

    if(mhd_glm) w(ixO^S,psi_)=0.d0

in subroutine `usr_init_one_grid` and ( subroutine `usr_special_bc` if exists).

#### Combined fix: typedivbfix='lindejanhunen' or 'lindepowel'

Combining diffusive fix and Janhunen or Powell fix by add both source
terms of these methods at the same time.
