# FUnTiDES — Acoustic Case

## 1. Variable Density Acoustic Wave Equation

Consider the variable density acoustics PDE parametrized in bulk modulus $\kappa$ and buoyancy $b = 1/\rho$.
This is an algebraically privileged parametrization, as can be seen in the derivation below:

$$\kappa \frac{\partial^2 p}{\partial t^2} - \nabla \cdot b \nabla p = f$$

where $p$ is the pressure field, $\kappa = \rho v_p^2$ is the bulk modulus, and $b = 1/\rho$ is the buoyancy.
The two model parameters targeted by inversion are:

- $\kappa$ — bulk modulus
- $b = 1/\rho$ — buoyancy

---

## 2. SEM Discretisation

The domain is decomposed into hexahedral elements. On each element $K^e$, the fields are expanded
on a tensor-product basis of Gauss–Lobatto–Legendre (GLL) polynomials of order $r$.
The GLL quadrature rule yields a diagonal mass matrix — the key property of SEM.

Each element has $(r+1)^3$ local quadrature/interpolation points, with local index
`lIdx = i + j*(r+1) + k*(r+1)²` mapping to a global node via `globalNodeIndex(e, i, j, k)`.

Given $p = \sum_{i=1}^{(r+1)^3} p_i \phi_i$, the semi-discrete (spatially discretised, continuous in time) system on element $K^e$ reads:

$$M^e(\kappa) \frac{p^{n+1} - 2p^n + p^{n-1}}{dt^2} + S^{\partial e}(\kappa,b) \frac{p^{n+1} - p^{n-1}}{2\,dt} + R^e(b)\, p^n = F^n$$

### 2.1 Element matrices

**Mass matrix** (diagonal for SEM):

$$M_{ij}(\kappa) = \int_{K^e} \kappa\, \phi_i \phi_j$$

**Stiffness matrix** (not diagonal):

$$R_{ij}(b) = \int_{K^e} b\, \nabla\phi_i \cdot \nabla\phi_j$$

**Boundary damping matrix** (on external boundary elements only):

$$S_{ij}\!\left(\sqrt{\tfrac{b}{\kappa}}\right) = \int_{\partial K^e_{\text{ext}}} \sqrt{\tfrac{b}{\kappa}}\, \frac{\partial \phi_i}{\partial n}\, \phi_j$$

Note that $M$ and $R$ are symmetric ($M_{ij} = M_{ji}$, $R_{ij} = R_{ji}$), while $S$ is not.
The following inner-product identities hold:

$$p^T M q = q^T M p, \qquad p^T R q = q^T R p, \qquad p^T S q = q^T S^T p$$

A key feature of this parametrization is that the matrices are linear in their parameter arguments:

$$\delta M_{ij}(\kappa) = M_{ij}(\delta\kappa), \qquad \delta R_{ij}(b) = R_{ij}(\delta b)$$

which makes the API for both the forward solver and gradient computation uniform — all kernels are
called with multiplier 1 for the gradient, and with the model value for the forward solve.

In the code, $R p$ is formed on-the-fly without storing $R$ and is referred to as the *stiffness vector*
(`INTEGRAL_TYPE::computeStiffnessTerm`). $M$ is assembled and stored as a diagonal vector
(`computeGlobalMassMatrix` / `massMatrixGlobal_`).

### 2.2 Leapfrog time integration

The solver uses a second-order leapfrog (Störmer–Verlet) scheme. At each global node $I$:

$$p^{n+1}_I = \frac{2 M_I(\kappa) p^n_I - \left(M_I(\kappa) - \tfrac{dt}{2} S_I\right) p^{n-1}_I - dt^2 F^n_I}{M_I(\kappa) + \tfrac{dt}{2} S_I}$$

where $S_I$ is the sponge/damping coefficient. Without damping this reduces to:

$$p^{n+1}_I = 2 p^n_I - p^{n-1}_I - \frac{dt^2}{M_I(\kappa)}\, F^n_I$$

The geometric nodal volume is:

$$\Omega_I = \sum_{e \ni I} w_I^e \, |J_I^e|$$

so that $M_I(\kappa) = \kappa_I \,\Omega_I$ (for node-based models) or $M_I(\kappa) = \kappa_e \,\Omega_I$ (element-based).

---

## 3. Forward Problem and Linearisation

Denote the forward operator by $L$, so $L p = f$. Its linearisation (Born approximation) is:

$$L\,\delta p = -\left( M^e(\delta\kappa) \frac{p^{n+1} - 2p^n + p^{n-1}}{dt^2} + R^e(\delta b)\, p^n \right)$$

which directly follows from the linearity $\delta M(\kappa) = M(\delta\kappa)$ and $\delta R(b) = R(\delta b)$.
The boundary term $\delta S$ is:

$$\delta S_{ij}\!\left(\sqrt{\tfrac{b}{\kappa}}\right) = S_{ij}\!\left( \frac{\delta b}{2\sqrt{\kappa b}} - \frac{\sqrt{b}\,\delta\kappa}{2\kappa\sqrt{\kappa}} \right)$$

---

## 4. Adjoint-State Gradient

### 4.1 Misfit functional

$$J(\kappa, b) = \frac{1}{2} \| \mathcal{R}\, p^n - d^n \|^2$$

where $\mathcal{R}$ is the receiver restriction operator and $d^n$ are the observed data.
The first variation is:

$$\delta J = \langle \delta p^n,\, \mathcal{R}^*(\mathcal{R} p^n - d^n) \rangle$$

### 4.2 Adjoint equation $L^* q = \mathcal{R}^*(\mathcal{R} p - d)$

The adjoint operator $L^*$ is derived from $\langle L p^n, q^n \rangle = \langle p^n, L^* q^n \rangle$
via summation by parts in time.

**Mass term** (symmetric, time-reversible):

$$\left\langle M^e \frac{p^{n+1} - 2p^n + p^{n-1}}{dt^2},\, q^n \right\rangle = \left\langle p^n,\, M^e \frac{q^{n-1} - 2q^n + q^{n+1}}{dt^2} \right\rangle \implies L_1^* = L_1$$

**Stiffness term** (symmetric):

$$\langle R\, p^n, q^n \rangle = \langle p^n, R\, q^n \rangle \implies L_2^* = L_2$$

**Boundary term** (not symmetric):

$$\langle S(p^{n+1} - p^{n-1}), q^n \rangle = \langle p^n, -S^T(q^{n-1} - q^{n+1}) \rangle \implies L_3^* \neq L_3$$

Boundary terms at initial/final time drop because $p^0 = p^1 = 0$ and $q^N = q^{N-1} = 0$.
Therefore:

$$L^* = L_1 + L_2 + L_3^*$$

The adjoint problem is a **final-value problem solved backwards in time**, with the same spatial
stiffness operator as the forward problem, but with a sign-flipped boundary damping term.

### 4.3 Gradient expressions

Introducing the adjoint field $q$ satisfying $L^* q = \mathcal{R}^*(\mathcal{R} p - d)$:

$$\delta J = \langle L\,\delta p,\, q \rangle = -\left\langle M^e(\delta\kappa) \frac{p^{n+1} - 2p^n + p^{n-1}}{dt^2} + R^e(\delta b)\, p^n,\; q_n \right\rangle$$

Expanding the buoyancy term explicitly and using symmetry of $R$:

$$\langle R(\delta b)\, p^n, q_n \rangle = \sum_n \sum_{ij} \delta b \int_{K^e} \nabla\phi_i \cdot \nabla\phi_j\, p_n^j\, q_n^i = \left\langle \delta b,\; \sum_n p_n^T R^e(1)\, q_n \right\rangle$$

The same approach applied to the mass term, using summation by parts to move the time
derivative from $p$ to $q$, gives:

$$\left\langle M^e(\delta\kappa) \frac{p^{n+1} - 2p^n + p^{n-1}}{dt^2}, q_n \right\rangle = \left\langle \delta\kappa,\; \sum_n p_n^T M^e(1) \frac{q^{n-1} - 2q^n + q^{n+1}}{dt^2} \right\rangle$$

The final gradient per element is:

$$\nabla J(\kappa, b) = -\begin{pmatrix} \displaystyle\sum_n \sum_{ij} \phi_i \phi_j \frac{q_{n-1}^i - 2q_n^i + q_{n+1}^i}{dt^2}\, p_n^j \\[8pt] \displaystyle\sum_n \sum_{ij} \nabla\phi_i \cdot \nabla\phi_j\; q_n^i\, p_n^j \end{pmatrix}$$

> **Practical note:** only $p$ needs to be saved during the forward pass. Both gradient terms
> are evaluated during the backward pass where $q$ is in memory.

### 4.4 Boundary contributions (from $S$)

For elements touching the external boundary, additional gradient terms arise from $\delta S$:

$$\nabla J(\kappa, b)\big|_{\partial K^e_{\text{ext}}} = -\begin{pmatrix} \displaystyle\sum_n p^T S^T\!\left(-\frac{\sqrt{b}}{2\kappa\sqrt{\kappa}}\right) \frac{q_{n-1} - q_{n+1}}{2\,dt} \\[8pt] \displaystyle\sum_n p^T S^T\!\left(\frac{1}{2\sqrt{\kappa b}}\right) \frac{q_{n-1} - q_{n+1}}{2\,dt} \end{pmatrix}$$

These are currently neglected in `DifferentiatorAcoustic` as they only affect boundary elements.

---

## 5. Discrete Gradient Assembly in Code

### 5.1 Model on elements (`IS_MODEL_ON_NODES = false`)

Each element has a unique gradient storage index — no race conditions.

**`grad_kappa[e]`** — mass-term integral over the element:

$$G^\kappa_e \mathrel{+}= \sum_q w_q |J^e_q|\; \frac{q_{n-1}^q - 2q_n^q + q_{n+1}^q}{dt^2}\; p_n^q$$

In code, $\ddot{q}$ is approximated on the fly as
$(\texttt{qnPrevPrev}[q] - 2\,\texttt{qnPrev}[q] + \texttt{qn}[q]) / dt^2$.

**`grad_buoyancy[e]`** — stiffness bilinear form over the element:

$$G^b_e \mathrel{+}= \sum_{i,j} R^e_{ij}(1)\; q_n^i\; p_n^j$$

No `ATOMICADD` needed; a thread-local accumulator writes once per element.

### 5.2 Model on nodes (`IS_MODEL_ON_NODES = true`)

Boundary/edge/corner nodes are shared between elements — `ATOMICADD` required.

**`grad_kappa[I]`** — scatter mass-term contribution to global node $I$:

$$G^\kappa_I \mathrel{+}= \sum_{e \ni I} w_I^e |J_I^e|\; \ddot{q}(I)\; p(I) = \Omega_I\; \ddot{q}(I)\; p(I)$$

**`grad_buoyancy[I]`** — scatter stiffness-term contribution to global node $I$ (test-function index):

$$G^b_I \mathrel{+}= \sum_{e \ni I} \sum_j R^e_{Ij}(1)\; q_n^j\; p_n^I$$

---

## 6. Preconditioning: From Raw Gradient to Smooth Sensitivity Kernel

### 6.1 What the accumulation produces

After summing over all time steps, the node-based raw gradient is:

$$G^\kappa_I = \Omega_I \cdot K^\kappa(x_I)$$

$\Omega_I$ varies strongly by topological position (corner vs. face vs. interior node), making
$G_I^\kappa$ **discontinuous** across element boundaries even when $K^\kappa$ is smooth.

### 6.2 Recovering the smooth kernel

The smooth $L^2$ sensitivity kernel is recovered by dividing by the geometric nodal volume:

$$K^\kappa(x_I) = \frac{G^\kappa_I}{\Omega_I}, \qquad K^b(x_I) = \frac{G^b_I}{\Omega_I}$$

$\Omega_I$ is assembled by running `computeMassTerm` with `model_factor = 1.0f` (no $\kappa$, no $\rho$).
It relates to the physical mass matrix as $\Omega_I = M_I(\kappa) / \kappa_I$.

### 6.3 Why the physical mass matrix $M_I(\kappa)$ must NOT be used

$$\frac{G^\kappa_I}{M_I(\kappa)} = \frac{\Omega_I \cdot K^\kappa(x_I)}{\Omega_I / \kappa_I} = \kappa_I \cdot K^\kappa(x_I)$$

The geometry cancels correctly, but the kernel is now multiplied by $\kappa_I$, artificially
amplifying the gradient where the model is large. The gradient direction in model space is wrong.

### 6.4 Where normalisation belongs

The normalisation $G_I / \Omega_I$ is time-independent. It should be applied **once** as a
post-processing step after the full time loop, at the optimisation/inversion layer:

```python
grad_kappa_smooth    = grad_kappa    / nodal_volume
grad_buoyancy_smooth = grad_buoyancy / nodal_volume
```

`DifferentiatorAcoustic` accumulates the raw gradient; preconditioning is the caller's responsibility.

### 6.5 Summary

| Quantity | Formula | Usage |
|---|---|---|
| Geometric nodal volume | $\Omega_I = \sum_{e \ni I} w_I^e \|J_I^e\|$ | Preconditioning denominator |
| Physical mass matrix | $M_I(\kappa) = \kappa_I\,\Omega_I$ | Leapfrog time update only |
| Raw `gradKappa[I]` | $\Omega_I \cdot K^\kappa(x_I)$ | Output of `DifferentiatorAcoustic` |
| Raw `gradBuoyancy[I]` | $\Omega_I \cdot K^b(x_I)$ | Output of `DifferentiatorAcoustic` |
| Smooth kernel | $K^\kappa(x_I) = G^\kappa_I / \Omega_I$ | Input to line search / model update |

---

## 7. Convolutional PML (C-PML) Absorbing Layer

The sponge layer of §2 is a *non-matched* absorber: it multiplies the field by a taper
$1/(1+\sigma)$ and reflects a fraction of the incident energy back into the domain. The
Convolutional PML (C-PML) replaces it with a *matched* absorbing layer that, for a plane
wave at normal incidence, is reflectionless at the interior interface and attenuates
exponentially inside the layer.

### 7.1 Stretched-coordinate formulation

In the PML layer the spatial derivative in direction $i$ is replaced by the stretched
derivative

$$s_i = \kappa_i + \frac{d_i}{\alpha_i + i\omega}, \qquad \frac{\partial}{\partial x_i} \to \frac{1}{s_i}\frac{\partial}{\partial x_i}$$

with $\kappa_i$ the coordinate-stretching factor, $d_i$ the damping profile and $\alpha_i$
the frequency-shift parameter. The second-order acoustic equation inside the layer
becomes

$$\frac{1}{\kappa}\frac{\partial^2 p}{\partial t^2} = \nabla \cdot \left(\frac{1}{\rho}\,\widetilde{\nabla p}\right)$$

where the stretched gradient is written in terms of a memory variable $\psi_i$
(Komatitsch & Martin 2007):

$$\widetilde{\nabla p}_i = \frac{1}{\kappa_i}\frac{\partial p}{\partial x_i} - \frac{1}{\kappa_i}\psi_i, \qquad
\frac{\partial \psi_i}{\partial t} + \left(\alpha_i + \frac{d_i}{\kappa_i}\right)\psi_i = \frac{d_i}{\kappa_i}\frac{\partial p}{\partial x_i}$$

The memory variables are advanced with the exact first-order convolution
(Wang, Lee & Teixeira 2006, eq. 21):

$$\psi_i^{n+1} = c^0_i\,\psi_i^n + c^1_i\left(\frac{\partial p}{\partial x_i}\right)^n, \qquad
c^0_i = e^{-(\alpha_i + d_i/\kappa_i)\,dt}, \qquad
c^1_i = \frac{d_i/\kappa_i}{\alpha_i + d_i/\kappa_i}\left(1 - c^0_i\right)$$

The weak form is obtained by multiplying by a test function $v$ and integrating by
parts. The implemented kernel uses the **two-sided** (weighted) form: both the trial
gradient and the divergence of the test function are stretched,

$$\int_\Omega \frac{1}{\kappa}\frac{\partial^2 p}{\partial t^2}\,v\,d\Omega
\;+\; \int_\Omega \left(\frac{1}{\kappa}\nabla v\right)\cdot\left(\frac{1}{\kappa}\,\widetilde{\nabla p}\right)d\Omega = 0$$

The divergence stretch introduces a second memory variable $\chi_i$ per direction,
advanced by the same convolution as $\psi_i$ but driven by the flux divergence
$\partial_i G_i$ (the reference divergence of the assembled flux $G = w\,\alpha\,\det
J\,J^{-1}\widetilde{\nabla p}$) instead of the pressure gradient. Both memory variables
are used at their *current* time level in the stretched quantities and advanced to the
next level afterwards, so the force at step $n$ uses $\psi^n,\chi^n$ — the standard
C-PML timing. Note that the $\frac{1}{\kappa}$ factor on the mass term is *not*
applied in the implementation (see §7.3); the stiffness operator is stretched on both
sides. This two-sided (weighted) form is a valid C-PML variant for the second-order
form — note that the canonical reference implementation (SPECFEM3D, Komatitsch et al.)
stretches only the trial gradient and leaves the divergence unstretched; both absorb,
and the two-sided form is what the sum-factorization kernel computes here.

### 7.2 Profiles

With $\delta$ the distance from the inner edge of the layer and $L$ its thickness, the
profiles are (Komatitsch & Martin 2007):

$$d_i(\delta) = d_{\max}\left(\frac{\delta}{L}\right)^N, \qquad
\kappa_i(\delta) = 1 + (\kappa_{\max}-1)\left(\frac{\delta}{L}\right)^N, \qquad
\alpha_i(\delta) = \alpha_{\max}\left(1 - \frac{\delta}{L}\right)^N$$

$$d_{\max} = -\frac{(N+1)\,v_p}{2L}\,\ln R$$

where $R$ is the target reflection coefficient, $N$ the profile exponent (default 2,
quadratic), $v_p$ the local P velocity, $\kappa_{\max}$ the maximum coordinate stretch
(default 1 = none) and $\alpha_{\max}$ the maximum frequency shift (default 0).

### 7.3 Implementation notes

- **Kernel.** `computeStiffnessTermSumFactPML` (both the makutu and tensorial backends)
  computes the reference gradient, the physical gradient $\nabla p = J^{-T}\nabla_\xi p$,
  builds the stretched gradient and assembles the flux
  $G = w\,\alpha\,\det J\,J^{-1}\widetilde{\nabla p}$ followed by the stretched
  divergence. The memory variables are advanced *inside* the stiffness kernel, so no
  separate memory-variable pass is needed.
- **Two-sided stretching.** Both the trial gradient and the divergence of the test
  function are stretched, each with its own memory variable:
  $$\widetilde{\nabla p}_i = \frac{1}{\kappa_i}\left(\frac{\partial p}{\partial x_i} - \psi_i\right), \qquad
  \widetilde{\mathrm{div}}_i = \frac{1}{\kappa_i}\left(\partial_i G_i - \chi_i\right)$$
  with $\psi_i$ (gradient memory) and $\chi_i$ (divergence memory) advanced by the same
  first-order convolution,
  $$\psi_i^{n+1} = c^0_i\,\psi_i^n + c^1_i\left(\frac{\partial p}{\partial x_i}\right)^n, \qquad
  \chi_i^{n+1} = c^0_i\,\chi_i^n + c^1_i\left(\partial_i G_i\right)^n$$
  Both are used at their *current* time level in the stretched quantities and advanced
  to the next level afterwards, so the force at step $n$ uses $\psi^n,\chi^n$ — the
  standard C-PML timing. With a zero profile ($d=0$, $\kappa=1$) $\psi_i$ and $\chi_i$
  stay zero and the stretched quantities equal the unstretched ones, so the PML kernel
  reduces exactly to `computeStiffnessTermSumFact`. This is the consistency test oracle.
- **Mass term.** The mass matrix is *not* stretched: the $\frac{1}{\kappa}$ factor on
  $\partial^2 p/\partial t^2$ in §7.1 is not applied in `computeMassTerm`. The
  implemented PML therefore stretches only the stiffness operator (both gradients),
  a valid C-PML variant; the mass stretch is a higher-order correction that is not
  needed for the absorbing behaviour.
- **Storage.** Coefficients are stored per node (`pmlCoefficients_`, 18 floats/node:
  $d,\kappa,\alpha,c^0,c^1,c^2$ per direction). The memory variables are stored per element
  per GLL point (`pmlMemoryVariables_`, $6\times(r+1)^3$ per element: $\psi$ and $\chi$,
  3 components each) because the gradient is element-local — a per-node storage would
  race on shared nodes.
- **Interaction with other boundaries.** In the PML region the sponge taper is set to 1
  (no double absorption) and the first-order absorbing BC of §2 is disabled (the PML
  replaces it). The free-surface condition is unaffected.
- **GEMM path.** The tensorial GEMM path precomputes $W = w\,\alpha\,B$, which is only
  valid for the unstretched operator. When the PML is enabled, the acoustic contribution
  is routed through the sum-factorization (Flat) kernel instead.
- **Adjoint mode.** The memory variables are advanced by the forward convolution inside
  the stiffness kernel. In backward/adjoint mode (`updateFieldsBackward`) the same kernel
  runs, but the convolution is not time-reversed, so the PML is not strictly
  adjoint-consistent there. Forward-mode PML absorption is unaffected.

### 7.4 CLI

```
--pml-size <meters>        PML thickness (0 = disabled, sponge used)
--pml-profile <N>          profile exponent (default 2)
--pml-reflection <R>       target reflection coefficient (default 1e-3)
--pml-alpha-max <a>        max frequency shift (default 0)
--pml-kappa-max <k>        max coordinate stretch (default 1)
```

The PML is only implemented for the acoustic physics; for elastic or acousto-elastic
simulations the flags are ignored.
