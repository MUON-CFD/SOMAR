The equations of motion
=======================

.. _eom_largescale:

The large-scale model
---------------------

To be verbose and specific, SOMAR solves the incompressible, Boussinesq
Navier-stokes equations in the following form.

.. math::
   :nowrap:

   \begin{align}
      \frac{\partial u}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i u
      +  2 J \left( \nu + \nu_{sgs} \right) S^{i x} \right]
      - \frac{\partial \xi}{\partial x} \frac{\partial p}{\partial \xi}
      + fv + F_u
      \\[1em]
      \frac{\partial v}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i v
      + 2 J \left( \nu + \nu_{sgs} \right) S^{i y} \right]
      - \frac{\partial \eta}{\partial y} \frac{\partial p}{\partial \eta}
      -fu + F_v
      \\[1em]
      \frac{\partial w}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i w
      + 2 J \left( \nu + \nu_{sgs} \right) S^{i z} \right]
      - \frac{\partial \zeta}{\partial z} \frac{\partial p}{\partial \zeta}
      - b' + F_w
      \\[1em]
      \frac{\partial T}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i T
      + J \left( \kappa_T + \kappa_{sgs,T} \right) \frac{\partial T'}{\partial \xi^i} \right]
      + F_T
      \\[1em]
      \frac{\partial S}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i S
      + J \left( \kappa_S + \kappa_{sgs,S} \right) \frac{\partial S'}{\partial \xi^i} \right]
      + F_S
      \\[1em]
      \frac{\partial \lambda}{\partial t} & =
        \frac{1}{J} \sum_i \frac{\partial}{\partial \xi^i}\left[ -U^i \lambda
      + J \left( \kappa_{\lambda} + \kappa_{sgs,\lambda} \right) \frac{\partial \lambda}{\partial \xi^i} \right]
      + F_{\lambda}
      \\[1em]
      0 & =
      \frac{\partial U}{\partial \xi}
      + \frac{\partial V}{\partial \eta}
      + \frac{\partial W}{\partial \zeta}
      \\[1em]
      b' & = \Gamma(T,S) - \Gamma(\overline{T}, \overline{S})
   \end{align}

The :math:`(\xi,\eta,\zeta)` are logical coordinates and the :math:`(x,y,z)` are physical coordinates. The logical and physical coordinates are equal unless using a stretched coordinate system. See `Coordinate systems`_ and :ref:`Stretching the grids<StretchingTheGrids>` for more details.

All primed quantities represent deviations from a static vertical
stratification. That is, :math:`T'(x,y,z,t) = T(x,y,z,t) - \overline{T}(z)`,
and so on. When a stratification is provided, splitting the active scalars
avoids diffusion of the background state and avoids repeated computation of the
dominant hydrostatic pressure component.

:math:`S^{i \alpha}` is the rate-of-strain tensor, made
curvilinear-based in the first index and Cartesian-based in the second.
Also, the capital velocity components represent the divergence-free
curvilinear-based advecting velocity components,
:math:`\left(U,V,W\right) = \left(J \frac{\partial \xi}{\partial x} u, J \frac{\partial \eta}{\partial y} v, J \frac{\partial \zeta}{\partial z} w \right)`.

By default, the custom forcing functions, :math:`F_{\bullet}` are zero, but these can be modified to fit your simulation's needs.
SOMAR also allows for custom passive scalars, :math:`\lambda`, which can be made active by modifying the custom forcing functions.
Although we only show one such scalar, SOMAR allows for any number of them.

.. warning::

    SOMAR uses a projection method to enforce incompressibility. This means the variable :math:`p` is technically not the physical pressure. It is whatever it needs to be to make the advecting velocity divergence-free. In the domain's interior, :math:`p` is a reasonable, low-order estimate of the pressure. However, it is computed using unphysical homogeneous Neumann boundary conditions, causing it to deviate from reality in boundary layers!


The subgrid-scale model
-----------------------

Without an LES, the sub-grid scale data :math:`\nu_{sgs},\,\kappa_{sgs,\bullet}`
are identically zero. When using an LES, these are computed as decribed in
`Ducros, et. al. <https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/largeeddy-simulation-of-transition-to-turbulence-in-a-boundary-layer-developing-spatially-over-a-flat-plate/C277DE968A1FD929D3CB05FDBC434AAD>`_

.. math::
   :nowrap:

   \begin{align}
     \nu_{sgs} &=
     0.0014\, C_K^{-3/2}\, \overline{\Delta}\, \sqrt{\overline{F_2}}
     \\[1em]
     \kappa_{sgs,\bullet} &= \nu_{sgs} / \text{Pr}_{sgs, \bullet}
     \\[1em]
     \overline{F_2} &=
     \frac{1}{6}
        \left[
            \left( \tilde{u}_{i+1,j,k} - \tilde{u}_{i,j,k} \right)^2
          + \left( \tilde{u}_{i-1,j,k} - \tilde{u}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{u}_{i,j+1,k} - \tilde{u}_{i,j,k} \right)^2
          + \left( \tilde{u}_{i,j-1,k} - \tilde{u}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{u}_{i,j,k+1} - \tilde{u}_{i,j,k} \right)^2
          + \left( \tilde{u}_{i,j,k-1} - \tilde{u}_{i,j,k} \right)^2
        \right] \\
     & +
     \frac{1}{6}
        \left[
            \left( \tilde{v}_{i+1,j,k} - \tilde{v}_{i,j,k} \right)^2
          + \left( \tilde{v}_{i-1,j,k} - \tilde{v}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{v}_{i,j+1,k} - \tilde{v}_{i,j,k} \right)^2
          + \left( \tilde{v}_{i,j-1,k} - \tilde{v}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{v}_{i,j,k+1} - \tilde{v}_{i,j,k} \right)^2
          + \left( \tilde{v}_{i,j,k-1} - \tilde{v}_{i,j,k} \right)^2
        \right] \\
     & +
     \frac{1}{6}
        \left[
            \left( \tilde{w}_{i+1,j,k} - \tilde{w}_{i,j,k} \right)^2
          + \left( \tilde{w}_{i-1,j,k} - \tilde{w}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{w}_{i,j+1,k} - \tilde{w}_{i,j,k} \right)^2
          + \left( \tilde{w}_{i,j-1,k} - \tilde{w}_{i,j,k} \right)^2 \right. \\ & \hspace{2em} \left.
          + \left( \tilde{w}_{i,j,k+1} - \tilde{w}_{i,j,k} \right)^2
          + \left( \tilde{w}_{i,j,k-1} - \tilde{w}_{i,j,k} \right)^2
        \right]

    \\[1em]
    \tilde{\varphi} &= \mathscr{L}^3\left[ \varphi \right]
    \\[1em]
    \mathscr{L}\left[\varphi\right] &=
      \left( \varphi_{i+1,j,k} -2 \varphi_{i,j,k} + \varphi_{i-1,j,k} \right) \\ &
    + \left( \varphi_{i,j+1,k} -2 \varphi_{i,j,k} + \varphi_{i,j-1,k} \right) \\ &
    + \left( \varphi_{i,j,k+1} -2 \varphi_{i,j,k} + \varphi_{i,j,k-1} \right),
   \end{align}

where :math:`C_K=0.5`,
:math:`\overline{\Delta}=(J \, \Delta\xi \, \Delta\eta \, \Delta \zeta)^{1/3}`,
and each scalar's subgrid Prandtl numbers are provided by the user through the input file.

When using multiple levels of adaptive refinement, the sub-grid parameterizations
:math:`\nu_{sgs}` and :math:`\kappa_{sgs,\bullet}` can be computed at all levels,
or only at the finest then interpolated down to coarser levels.


Coordinate systems
------------------

.. _eom_coordsys:

SOMAR bounces back and forth between two coordinate systems, *logical* and *physical*. The logical coordinates of each node are always defined by

.. math::

    \xi_i   &= i \Delta\xi   &= i L_x/N_x, \\
    \eta_j  &= j \Delta\eta  &= j L_y/N_y, \\
    \zeta_k &= k \Delta\zeta &= k L_z/N_z.

Where :math:`(i,j,k)` are the nodal indices, :math:`L_{\bullet}` is the domain's length and :math:`N_{\bullet}` is the number of cells in each direction. The logical coordinates at cell-centers are defined similarly, but with 1/2 added to the indices.

The physical coordinates are computed from the logical coordinates via three stretching functions,

.. math::

    \left( x, y, z \right) = \left( f(\xi), g(\eta), h(\zeta) \right).

The identity map, :math:`\left( x, y, z \right) = \left( \xi, \eta, \zeta \right)`, leads to a simple, unstretched Cartesian coordinate system. The user can customize these stretching functions to meet the needs of the simulation.

The Jacobian matrix is :math:`\rm{diag}\left( \frac{\partial x}{\partial \xi}, \frac{\partial y}{\partial \eta}, \frac{\partial z}{\partial \zeta} \right)` and it's determinant is :math:`J = \frac{\partial x}{\partial \xi} \frac{\partial y}{\partial \eta} \frac{\partial z}{\partial \zeta}`. With this system, volume integrals of scalar fields are, for example,

.. math::

    M = \iiint \rho\,J\,d\xi\,d\eta\,d\zeta = \sum_{i,j,k} \left( \rho \, J \right)_{i,j,k} \, \Delta\xi \Delta \eta \Delta \zeta.

.. Spatial derivatives of scalar fields are, for example,

.. .. math::

..     \frac{\partial p}{\partial x} = \frac{\partial \xi}{\partial x} \frac{\partial p}{\partial \xi}.

.. In SOMAR, we use the word *flux* to mean a vector that has been scaled by :math:`J`. For example, the advective flux of a scalar, :math:`s`, is not :math:`\vec{u} s` but :math:`J \vec{u} s := \vec{U}s`. Also, we define the operators :math:`\mathbf{G}[\,] := J\,\nabla[\,]` and :math:`\mathbf{D}[\,] := J^{-1}\,\nabla \cdot [\,]`. Despite the differences in scaling, we simply call these the gradient and divergence operators. Written in full, the gradient of a scalar, :math:`s`, is

.. .. math::

..     \mathbf{G}s := J\,\nabla s = J \left( \frac{\partial \xi}{\partial x} \frac{\partial s}{\partial \xi}, \, \frac{\partial \eta}{\partial y} \frac{\partial s}{\partial \eta}, \, \frac{\partial \zeta}{\partial z} \frac{\partial s}{\partial \zeta} \right),

.. and the divergence of a flux-valued vector field, :math:`\vec{U} = (U,V,W)`, is

.. .. math::

..     \mathbf{D}\vec{U} := \frac{1}{J} \nabla \cdot \vec{U} = \frac{1}{J} \left( \frac{U}{\ \right)
