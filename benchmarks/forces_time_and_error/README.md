# Benchmark: Forces time and error

Here I'll do some benchmarks of the octree, looking both to the error (compared to the direct forces) and the time to compute the forces. I'll also consider the three implementations of multipoles: monopole (1), quadrupole (4) and octupole (8).

I also reproduce here some results of McMillan and Aarseth (1993), where they implement the octupole for collisional stellar systems.

---

## McMillan and Aarseth Results

As expected.

<img src="img/aarseth/aarseth_result.png" width="400px">

---

## Multipoles

<img src="img/multipoles_error_x_time/multipoles.png" width="400">

---

## Bigger $N$

<img src="img/bigger_N/times_BIG_Monopole_nt25.png" width="400">
<img src="img/bigger_N/times_BIG_Quadrupole_nt25.png" width="400">
<img src="img/bigger_N/times_BIG_Octupole_nt25.png" width="400">

---

# References

**[1]** McMillan, S. L. W. and Aarseth, S. J., “An O(N N) Integration Scheme for Collisional Stellar Systems”, <i>The Astrophysical Journal</i>, vol. 414, IOP, p. 200, 1993. doi:10.1086/173068.