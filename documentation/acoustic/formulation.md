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
