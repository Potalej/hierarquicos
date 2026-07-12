# How to generate initial values

To generate initial values like another N-body-ers usually do, we can use the Galactic Python [*galpy*](https://docs.galpy.org/en/stable) package.

## Plummer density profile

The most common density profile is the *Plummer profile* (Plummer, 1911) given by

$$
\rho(r) = \dfrac{3 M b^2}{4 \pi} (r^2 + b^2)^{5/2},
$$

where $M$ is the total masses of the system and $b$ is the *Plummer length/scale*. It relates to the half-mass radius $r_h$ (i.e. the radius that contains a total mass of $M/2$) as $r_h \approx 1.3 b$.

The correspondent potential is given by:

$$
\Phi(r) = - \dfrac{G M}{\sqrt{r^2 + b^2}},
$$

where $M$ is the total mass of the system. The use of a $\varepsilon$ as a softening to the Newton's classical potential is related directly to this density profile.

To generate random radius $r$ with this density we can just define the CDF (in spherical coordinates):

$$
P(r) = \dfrac{M(\lt r)}{M} = \dfrac{4 \pi}{M} \int_0^r \rho(q) q^2 dq,
$$

generate a random $X_1 \sim U[0,1]$ and solve $P(r) = X_1$. In this case specifically we can obtain $P^{-1}(X_1)$ explicitly:

$$
r = P^{-1}(X_1) = \dfrac{b}{\sqrt{X_1^{-2/3} - 1}}.
$$

> (Aarseth, 2003) suggests to use a rejection parameter to reject rare large values (e.g. $r > 10 r_h$).


We can guarantee that we are generating values with the correct distribution we can obtain the PDF for the Plummer profile:

$$
p(r) = \dfrac{dP(r)}{dr} = \dfrac{3 b^2 r^2}{(b^2 + r^2)^{5/2}}.
$$

To generate the positions $(x,y,z)$ we just generate two more $X_2, X_3 \sim U[0,1]$ and define

$$
z = (1 - 2 X_2)r
$$

$$
x = \sqrt{r^2 - z^2} \cos{2 \pi X_3},
$$

$$
y = \sqrt{r^2 - z^2} \sin{2 \pi X_3}.
$$

## Velocities and some very funny theory

To generate velocities we need to choose the way we want the system to evolve. We'll think that's ok that no direction is the favorite for the velocity, so the system is *isotropic*. This is apparently the most common and the simpler to.

### Collisionless Boltzmann equation

First, let's consider a distribution function $f(\vec x, \vec v)$ for the bodies that relates to $\rho$ as:

$$
\rho(\vec x) = \int d^3 v \ f(\vec x, \vec v).
$$

As long as no body should just appear ou disappear, $f$ represents a conservative field, so it's flow preserves volume. By the Liouville Theorem we get

$$
\dfrac{d f}{dt} = \dfrac{\partial f}{\partial t} + \dfrac{d \vec x}{dt} \dfrac{\partial f}{\partial \vec x} + \dfrac{d \vec v}{dt} \dfrac{\partial f}{\partial \vec v} = 0.
$$

This is the *Continuity Equation*, or as it's know in the streets: *the collisionless Boltzmann equation*.

In fact, as long as we are in a Hamiltonian system and all this thing can be done using $(\vec q, \vec p)$ instead of $(\vec x, \vec v)$, this is just the equation for an observable under the Hamiltonian system using the Poisson bracket:

$$
\dfrac{\partial f}{\partial t} + \{H, f\}_{(\vec p, \vec q)} = 0.
$$

### Jeans Theorem

It's important to note that if we consider a system under equilibrium (in some sense), $\partial f / \partial t = 0$, therefore what we have is the Jeans Theorem: *any function of the integrals of motion is a solution of the collisionless Boltzmann equation in equilibrium, and any solution of this equation only depends on the integrals of motion*. The proof is trivial at this point.

### Ergodic distribution functions

We this Theorem we can define the *relative potential* and the *relative/specific energy*, respectively, as:

$$
\Psi = - \Phi + \Phi_0,
\quad
\mathcal E = - E + \Phi_0 = \Psi - \dfrac{1}{2} v^2,
$$

where $\Phi_0$ is some convenient constant and $v = |\vec v|$. It would be good to find some function $f$ that only depends on $\mathcal E$, as they have a cool name: *ergodic distribution function*. In this case we have:

$$
\rho(r) = \int d \vec v \ f(r, \vec v) = 4 \pi \int dv \ v^2 f(\mathcal E)
= 4 \pi \int_0^\Psi d \mathcal E \ \sqrt{2 (\Psi - \mathcal E)} f(\mathcal E).
$$

The potential is spheric and therefore monotonic radially, so we can write $\rho$ as a function of $\Psi$ in a well-defined way:

$$
\dfrac{1}{\sqrt{8} \pi} \rho(\Psi) = 2 \int_0^\Psi d \mathcal E \ \sqrt{\Psi - \mathcal E} \ f (\mathcal E).
$$

We then differentiate it wrt to $\Psi$:

$$
\dfrac{1}{\sqrt{8} \pi} \dfrac{d \rho(\Psi)}{d \Psi} = \int_0^\Psi d \mathcal E \dfrac{f(\mathcal E)}{\sqrt{\Psi - \mathcal E}}.
$$

This is an Abel integral equation (btw I don't have any idea of what it means), and it can be inverted to

$$
f(\mathcal E) = \dfrac{1}{\sqrt{8} \pi^2} \dfrac{d}{d \mathcal E} \int_0^\mathcal{E} d \Psi \ \dfrac{1}{\sqrt{\mathcal E - \Psi}} \dfrac{d\rho}{d \Psi}.
$$

This is the well-called **Eddington inversion formula**. Is not trivial to evaluate this thing for any density profile, but for the Plummer it is. We get:

$$
f(\vec x, \vec v) = \dfrac{24 \sqrt{2}}{7 \pi^3} \dfrac{b^2}{G^5 M^4} (-E(\vec x, \vec v))^{7/2} \equiv f(E).
$$

Considering a volume element in the velocity space between $v$ and $dv$, we have 

$$
p(v | r) \propto 4 \pi v^2 f(E) dv.
$$

Using the value of $E$, we get

$$
p(v | r) \propto v^2 \left(-\dfrac{v^2}{2} - \Phi(r)\right)^{7/2}.
$$

Let's suppose that we don't want particles with a velocity bigger than the escape velocity $v_{esc}(r) = \sqrt{- 2 \Phi(r)}$. Defining $q := v/v_{esc} \in [0,1]$, we can rewrite the probability:

$$
p(v | r) \propto q^2 (1 - q^2)^{7/2} =: g(q).
$$

To generate a velocity given $r$, we need to use rejection (von Neumann). First, $g(q)$ has a global maximum (as $q \geq 0$) at $q_{max} = \sqrt{2}/3$, with $g(q_{max}) \approx 0.0922 $. Now we do:

1. Generate $q \sim U[0,1]$.
2. Generate $y = g(q_{max}) \cdot y_0$, with $y_0 \sim U[0,1]$.
3. If $y \leq g(q)$, we accept and define $v = q \ v_{esc}$. If not, go to (1).

> Very nice!

### Spherical Jeans equation and the equilibrium

Anyway, we can do some *maracutaia* to write the collisionless Boltzmann equation in spherical coordinates, multiply it by the radial momentum, integrate it over the 3-dimensional momentum and vanishes some terms basing on the Jeans Theorem. What we got is the **Spherical Jeans Equation**

$$
\dfrac{d(\rho \sigma_r^2)}{dr} + \dfrac{2 \beta \rho \sigma_r^2}{r} = - \rho \dfrac{d \Phi}{dr},
$$

where $\sigma_r$ is the mean radial velocity and $\beta$ is the orbital anisotropy measure, given by:

$$
\beta = 1 - \dfrac{\sigma_\theta^2 + \sigma_\varphi^2}{2 \sigma_r^2}.
$$

Note: if a system is isotropic, $\sigma_r = \sigma_\theta = \sigma_\varphi$, therefore $\beta = 0$ and the Jeans Equation reduces to an ODE:

$$
\dfrac{d (\rho \sigma_r^2)}{dr} = - \rho \dfrac{d \Phi}{dr}.
$$

Again, using the monotonicit of the potential $\Phi$, it's reasonable to suppose $\Phi(\infty) = 0$. With this, we can multiply both sides of the ODE by $4 \pi r^3$ (*why not?*) and integrate over $r$:

$$
\int_0^\infty 4 \pi r^3 \dfrac{d (\nu \sigma_r^2)}{dr} dr = - \int_0^\infty 4 \pi r^3 \rho \dfrac{d \Phi}{dr} dr.
$$

With integration by parts and some definitions, we can get for the left side:

$$
\int_0^\infty 4 \pi r^3 \dfrac{d (\nu \sigma_r^2)}{dr} dr = - 12 \pi \int_0^\infty r^2 \rho(r) \sigma_r^2 dr = - 2 T,
$$

the kinect energy!

For the right side we use that:

$$
\dfrac{d\Phi}{dr} = \dfrac{G M(\lt r)}{r^2}.
$$

We get:

$$
-4 \pi \int_0^\infty r^3   \rho(r) \dfrac{d \Phi}{dr} dr = - 4 \pi G \int_0^\infty r \rho (r) M(\lt r) dr,
$$

but $M(\lt r) = - r \Phi(r) / G$, then:

$$
= 4 \pi \int_0^\infty r^2 \rho(r) \Phi(r) dr
= \int \rho(\vec x) \Phi(\vec x) d^3 \vec x,
$$

and this is exactly the potential energy $V$ of the system!

Therefore, what we got is just the *virial theorem*: $2 T = -V$. This means that when we generate initial values using all this thing and suppose isotropic velocities, these values are in initial equilibrium (or at least near the equilibrium).

For the Plummer profile, we can evaluate either explicitly:

$$
V = \dfrac{1}{2} \int \rho(\vec x) \Phi(\vec x) d^3x = 2 \pi \int_0^\infty \rho(r) \Phi(r) r^2 dr = - \dfrac{3 \pi}{32} \dfrac{G M^2}{b}.
$$

$$
T = \dfrac{3 \pi}{64}\dfrac{G M^2}{b}
$$

and the total energy is:

$$
E = - \dfrac{3 \pi}{64}\dfrac{G M^2}{b}.
$$

Obviously this is exact only with $N \to \infty$, but we can expected something nearly these values with finite $N$.


## Using galpy

To sample positions following the Plummer density with Galpy we need the modules `potential` and `df`:

```python
from galpy import potential, df

def plummer (N, b=1.0):
  pot = potential.PlummerPotential(amp=1.0, b=b, normalize=False)
  dfunc = df.isotropicPlummerdf(pot=pot)
  orbits = dfunc.sample(n=N)
  return orbits
```

With this variable `orbits` we can extract the positions and velocities by
```python
xs,   ys,  zs = orbits.x(),  orbits.y(),  orbits.z()
vxs, vys, vzs = orbits.vx(), orbits.vy(), orbits.vz()
```

![](./plummer_radius_example_without_cutoff.png)

> In fact, the cut-off parameter may be necessary to avoid outliers and get something bounded. Applying a rejection parameter we get "better" numbers:

![](./plummer_radius_example_cutoff.png)

To verify the velocities we need to get the random variables $q = v / v_{esc}$ instead of just $v$.

```python
rs = np.linalg.norm(qs, axis=1)
pot = potential.PlummerPotential(amp=1.0, b=b, normalize=False)
v_esc = np.sqrt(-2.0 * pot(rs,0))
qs = vs / v_esc
```

Normalizing the values to compare with the PDF, we get:

![](./plummer_velocity_example.png)

For the dynamic properties, we can see that $N=10^3$ is not really so much, but we have something near the expected:

```
Dynamical properties (||real - expected||):
Potential V: 0.006832137175154629
Kinect T:    0.00337158496224943
Total E:     0.010203722137404059
Virial Q:    0.022781918828791548
```

---

## References 

* Aarseth Sverre J. Gravitational N-Body simulations. Cambridge University Press (2003)
* Bovy J. Dynamics and Astrophysics of Galaxies. Princenton University Press (2026)
* Plummer H. On the problem of distribution in globular star clusters. Mon. Not. Roy. Astron. Soc. 71, 460–470 (1911).
