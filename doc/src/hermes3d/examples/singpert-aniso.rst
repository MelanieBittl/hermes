Singular Perturbation
=====================

**Git reference:** Examples `singular perturbation <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes3d/examples/singpert-aniso>`_.

We solve a singularly perturbed elliptic problem that exibits a thin anis-tropic boundary layer 
that is difficult to solve. This examples demonstrates how the anisotropic refinements can save 
a big amount of degrees of freedom.

.. index::
   single: mesh; dynamical
   single: problem; elliptic

The computational domain is the unit cube $(0, 1)^3$, and the equation solved has the form :

.. math::
   :nowrap:
   :label: sing-perturb

   \begin{eqnarray*}
   - \Delta u + K^2 u &= K^2 &\hbox{ in }\Omega \\
                    u &= 0 &\hbox{ on }\partial\Omega
   \end{eqnarray*}

The boundary conditions are homegeneous Dirichlet. The right-hand side is chosen to keep the 
solution $u(x,y,z) \approx 1$ inside the domain. Here, we choose $K^2 = 10^4$ but everything 
works for larger values of $K$ as well. 

It is quite important to perform the initial refinements towards the boundary, thus providing 
a better initial mesh for adaptivity and  making convergence faster. 

Convergence graphs:

.. image:: singpert-aniso/singpert-aniso-conv.png
   :scale: 50%


.. image:: singpert-aniso/singpert-aniso-conv-time.png
   :scale: 50%



Solution and hp-mesh:

.. image:: singpert-aniso/singpert-aniso-sln.png
   :scale: 50%


.. image:: singpert-aniso/singpert-aniso-order.png
   :scale: 50%


