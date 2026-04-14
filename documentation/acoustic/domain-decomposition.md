# SEM MPI exchange options

## Governing Equation and Splitting

We start from the semi-discrete central-difference scheme.

The mass matrices are **diagonal**, so their inverses are applied **component-wise** (a scalar per degree of freedom).


$$
\quad M\,(p_{n+1}-2p_n+p_{n-1}) = \Delta t^{2}\,(Kp_n + s)
$$
$$
\Rightarrow \quad p_{n+1} = 2p_n - p_{n-1} + \Delta t^2\, M^{-1} (K p_n + s)
$$

## MPI border case

Then we split into subdomains and define:
- $r$ a subdomain MPI rank
- $M_r$ the element local mass matrix computed for rank $r$
- $K_r$ the element local stiffness matrix for rank $r$
- $n$ the number of neighbouring subdomains contributions for the mass matrix around the element
- $\Gamma$ the boundary between subdomains

$$
M = \sum_{r=0}^{n} M_r,\; K = \sum_{r=0}^{n} K_r.
$$

$$
M_r = 0, \quad K_r = 0 \quad \text{if element} \notin \Gamma
$$


We define $R_r$ the contribution of the rank $r$ to the global solution on the border of the subdomain.

$$
R_r = 2p_n - p_{n-1} + \Delta t^2 ( \frac{\sum M_r}{r})^{-1} (K_r p_n + s),
$$
$$
p_{n+1} = \frac{\sum R_r}{r}
$$

## 1D Example

Consider a 1D domain split between two MPI ranks:

```
Rank 0:           Rank 1:
[0---1---2---3] | [3---4---5---6]
        ^
        Γ (shared node 3)
```

At the shared node 3:
- Rank 0 computes: `M_0[3]`, `K_0[3]` (from elements [2-3])
- Rank 1 computes: `M_1[3]`, `K_1[3]` (from elements [3-4])

Local contributions:
$$
R_0 = 2p_n[3] - p_{n-1}[3] + \Delta t^2 (\frac{M_0[3]+M_1[3]}{2})^{-1} (K_0[3] p_n + s)
$$
$$
R_1 = 2p_n[3] - p_{n-1}[3] + \Delta t^2 (\frac{M_0[3]+M_1[3]}{2})^{-1} (K_1[3] p_n + s)
$$

Global update at node 3:
$$
p_{n+1}[3] = \frac{R_0 + R_1}{2}
$$

## 2D Example

Consider a 2D domain split between four MPI ranks in a 2×2 configuration:

```
┌─────────────────┬─────────────────┐
│                 │                 │
│    Rank 0       │     Rank 1      │
│                 │                 │
│                 │                 │
├─────────────────●─────────────────┤  ← Γ (horizontal boundary)
│                 │                 │
│    Rank 2       │     Rank 3      │
│                 │                 │
│                 │                 │
└─────────────────┴─────────────────┘
                  ↑
                  Γ (vertical boundary)
                  ● = shared corner node
```

At the shared corner node ● (center point where all 4 subdomains meet):
- Rank 0 computes: `M_0[●]`, `K_0[●]` (from elements in the bottom-right corner of Rank 0)
- Rank 1 computes: `M_1[●]`, `K_1[●]` (from elements in the bottom-left corner of Rank 1)
- Rank 2 computes: `M_2[●]`, `K_2[●]` (from elements in the top-right corner of Rank 2)
- Rank 3 computes: `M_3[●]`, `K_3[●]` (from elements in the top-left corner of Rank 3)

Local contributions:
$$
R_0 = 2p_n[●] - p_{n-1}[●] + \Delta t^2 \left(\frac{M_0[●]+M_1[●]+M_2[●]+M_3[●]}{4}\right)^{-1} (K_0[●] p_n + s)
$$
$$
R_1 = 2p_n[●] - p_{n-1}[●] + \Delta t^2 \left(\frac{M_0[●]+M_1[●]+M_2[●]+M_3[●]}{4}\right)^{-1} (K_1[●] p_n + s)
$$
$$
R_2 = 2p_n[●] - p_{n-1}[●] + \Delta t^2 \left(\frac{M_0[●]+M_1[●]+M_2[●]+M_3[●]}{4}\right)^{-1} (K_2[●] p_n + s)
$$
$$
R_3 = 2p_n[●] - p_{n-1}[●] + \Delta t^2 \left(\frac{M_0[●]+M_1[●]+M_2[●]+M_3[●]}{4}\right)^{-1} (K_3[●] p_n + s)
$$

Global update at node ●:
$$
p_{n+1}[●] = \frac{R_0 + R_1 + R_2 + R_3}{4}
$$

Note that:
- Edge nodes (e.g., node 17 or 10) are shared between 2 ranks
- Corner nodes (e.g., node 24) are shared between 4 ranks
- The number of ranks `n_ranks` in the averaging step depends on the node's position at the boundary

## Development

### Pseudo-Code

$$
\begin{aligned}
&\rule{10cm}{0.4pt} \\
&\textbf{Initialization phase (before time-stepping)} \\
\\
&\quad \text{! Exchange and average mass matrices at boundaries (pysabl).} \\
&\quad \textbf{for each } i \in \Gamma_r: \\
&\quad\quad M_r[i] \leftarrow \text{MPI\_Allreduce}(M_r[i], \text{MPI\_SUM}, \text{neighbors}(i)) \\
&\quad\quad n_{\text{ranks}} \leftarrow \text{count\_sharing\_ranks}(i) \\
&\quad\quad M_r[i] \leftarrow M_r[i] / n_{\text{ranks}} \\
\\
&\quad \text{! Compute inverse mass matrix for all nodes (funtides).} \\
&\quad M_r^{-1} \leftarrow 1.0 / M_r \\
&\rule{10cm}{0.4pt} \\
&\textbf{Time-stepping loop} \\
\\
&\quad \textbf{for } n = 0 \text{ to } N_{\text{steps}}: \\
&\quad\quad \text{! Compute local wavefield contribution } R_r \text{ for all nodes (funtides)} \\
&\quad\quad p_{n-1} \leftarrow 2p_n - p_{n-1} + \Delta t^2 \, M_r^{-1} \, (K_r \, p_n + s) \\
\\
&\quad\quad \text{! Swap time levels (fundtides)} \\
&\quad\quad p_{n} \leftarrow p_{n-1} \\
&\quad\quad p_{n-1} \leftarrow p_{n} \\
\\
&\quad\quad \text{! Exchange and average boundary contributions (pysabl)} \\
&\quad\quad \textbf{for each } i \in \Gamma_r: \\
&\quad\quad\quad p_{n}[i] \leftarrow \text{MPI\_Allreduce}(p_{n}[i], \text{MPI\_SUM}, \text{neighbors}(i)) \\
&\quad\quad\quad n_{\text{ranks}} \leftarrow \text{count\_sharing\_ranks}(i) \\
&\quad\quad\quad p_{n}[i] \leftarrow p_{n}[i] / n_{\text{ranks}} \\
&\rule{10cm}{0.4pt} \\
\end{aligned}
$$

### Complexity Analysis

| Algorithm Step | Memory Complexity | Communication Complexity |
|---------------|------------------|-------------------------|
| Mass matrix exchange (initialization) | No additional buffers | $O(\|\Gamma_r\| \cdot \log(n_{\text{neighbors}}))$ |
| Wavefield exchange (per time step) | No additional buffers | $O(\|\Gamma_r\| \cdot \log(n_{\text{neighbors}}))$ |

Where  $|\Gamma_r|$ is the boundary size for the rank $r$.

### Required changes

funtides:
- expose the mass matrix $M_r$ via pybind
- decouple the inversion of mass matrix from the solver intialisation
- expose the mass matrix inversion via pybind

sabl:
- no changes required

makutu-wave:
- reduce average the mass matrix at boundaries
- reduce average the wavefield at boundaries