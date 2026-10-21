# knotgridrand

Python tools for generating, representing, randomising and analysing knots on the cubic lattice.

## Overview

`knotgridrand` is a Python project developed for computational experiments with lattice knots. It provides tools for working with knots represented as ordered sequences of three-dimensional lattice coordinates and generating randomised knot configurations.

The package was developed as part of an MSci project investigating [**Stochastic Randomisation And Neural Network
Based Classification Of Knots In 3D Space**](https://drive.google.com/file/d/10OT540xMB-NDDF0PgPgucR_COFQMPtku/view?usp=sharing).

## Features

* Represent knots as ordered 3D lattice coordinates
* Convert between coordinate arrays and dictionary-based representations
* Generate collections of base knot configurations
* Randomise knots using type-preserving moves
* Run the BFACF algorithm[^bfacf] in three dimensions
* Introduce random large scale deformations using the Slide move (which translates a strand)
* Introduce large scale bends using consecutive Slide moves
* Calculate properties of knot confirgurations (length, energy, COM)
* Calculate and work with knot writhe
* Generate large datasets of randomised knots for machine-learning applications

[^bfacf]: EJ Janse Van Rensburg and SG Whittington. “The BFACF algorithm and knot-
ted polygons”. In: Journal of Physics A: Mathematical and General 24.23 (1991),
p. 5553.
## Knot Representation

Knots are represented as ordered points on a three-dimensional cubic lattice:

```text
(x₀, y₀, z₀)
(x₁, y₁, z₁)
...
(xₙ, yₙ, zₙ)
```

Consecutive points are connected by unit lattice steps along the Cartesian axes. Diagonal steps are not permitted.

A dictionary-based representation is used to make it easier to perform local operations and efficiently identify neighbouring lattice points.

The dictionary representation stores the knot points as {index : x,y,z}, with forward and backward connectivity maps and a {'x,y,z' : index} dictionary used for collision detection.

## Randomisation

The randomisation process relies on the **BFACF algorithm** for local Monte-Carlo random walks, the Slide move to introduce large scale deformations and a 'potential energy' that is kept track of to control the probability of accepting a move.

An energy function can be used to control the randomisation:

$$
E=\alpha N + \beta \left(N - \sum_{n=1}^{N}{a_n \cdot b_n} \right) + \gamma \left(\sum_{n=1}^{N}{r_n - r_{CoM}} \right)
$$

where:
* $N$ is the number of segments/points of the knot
* $a_n$ and $b_n$ are the two vectors connecting the point $n$ to the preceding and proceeding points
* $r_n$ and $r_{CoM}$ are the position vectors of point $n$ on the knot and the centre of mass, respectively 
* $\alpha$, $\beta$ and $\gamma$ control the relative energetic costs of knot length, number of bends and distance of each point to the CoM, respectively

A proposed move $M$ is accepted with probability $P(M)$:

$$
P(M) = 
\begin{cases}
e^{-\frac{\Delta E}{T}}, & \Delta E > 0, \\
1 - e^{-\frac{\Delta E}{T}}, & \Delta E \leq 0
\end{cases}
$$

where:
* $\Delta E$ is the change in energy between the two configurations
* $T$ is a temperature-like parameter 
## Writhe

The writhe of a knot can be calculated from its three-dimensional coordinates according to:

$$
    \omega_{StA}(i) = \frac{2}{4\pi}\sum_{k=i-\frac{l_w}{2}}^{k=i-\frac{l_w}{2}}
    \sum_{l=0}^{N}{\frac{(\mathbf{t}_k\times \mathbf{t}_l)\cdot
    (\mathbf{r}_k-\mathbf{r}_l)}{|\mathbf{r}_k-\mathbf{r}_l|^3}}
$$

[^writhe]
[^writhe]: From the Supplementary Information for Joseph Lahoud Sleiman et al. “Geometric learning of knot topology”. In: Soft Matter
20.1 (2024), pp. 71–78.

where:
* $t_k$ is the directional vector between the points $k$ and $k+1$ and is tangential to that segment of the knot (the same goes for $t_l$)  
* $\mathbf{r}_k-\mathbf{r}_l$ is the vector connecting every pair of points $k$ and $l$ along the curve

## Installation

Clone the repository:

```bash
git clone https://github.com/Ahmedg2v/knotgridrand.git
cd knotgridrand
```

(Optional) Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Windows:

```bash
.venv\Scripts\activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage

[examples.py](examples.py) contains examples of how the different functionalities in `knotgridrand` might be used.

More examples will be added as time passes/developments are made.

## Scientific Applications

The package was originally developed to support the generation of training data for a machine-learning project on knot classification.

In the project, large numbers of randomised knot configurations were generated computationally and subsequently used as input to neural-network models.

The workflow involved:

```text
Base knots
    ↓
Lattice representation
    ↓
Randomisation process
    ↓
Randomised knot configurations
    ↓
Geometric/topological features or raw coordinates
    ↓
Machine-learning classification
```

The project is therefore intended primarily as a computational research tool, but its individual components can also be used for experiments involving lattice knots, Monte Carlo methods, and computational knot theory.

## Performance

The package has been used to generate datasets containing **hundreds of thousands of randomised knot configurations**, with individual configurations containing approximately **1,000 lattice points on average**, using the University of Birmingham HPC cluster (BlueBear).

Each dataset of a **hundred thousand randomised knots** (of the same type) was generated using the scripts in QQQQQQQQ as an array job with 1000 nodes, randomising 100 knots each. This process took 15-20 minutes per dataset using the parameters referenced in the project REPORT.

## Dependencies

The project was developed in Python 3.10 and primarily uses the scientific Python ecosystem, including:

* NumPy
* SciPy
* Matplotlib

See `requirements.txt` for the complete list of dependencies.

## Acknowledgements

This project was completed as part of the MSci Theoretical Physics programme at the University of Birmingham under the supervision of Dr Jonathan S. Watkins.

## Author

**Ahmed Ismaiel**

## License

This project is currently available for personal and academic use. See the repository for licensing information.


