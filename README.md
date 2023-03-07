# PaScaL_TCS

Parallel and Scalable Library for Turbulent Convection Solver

The PaScal TCS is an efficient and scalable solver for large-scale direct numerical simulations of turbulent convective
flows, considering temperature-dependent fluid properties. In order to increase scalability on a massive-scale
distributed memory system, the PaScal_TCS decomposes the computational domain into a cubic-subdomain. The numerical
procedure is based on the monolithic projection method with staggered time discretization (MPM-STD), which is 
an implicit and non-iterative solver for wall-bounded turbulent flows. 
The PaScaL_TCS has the following features.

 1. Temperature-dependent fluid properties considering non-Oberbeck–Boussinesq effect
 2. Three-dimensional domain decompositions to exploit more parallelism
 3. PaScaL_TDMA library to solve the batched tri-diagonal systems partitioned by the domain decomposition. 
 4. Two transpose schemes in parallel fast Fourier transform (FFT) for the pressure Poisson equation solver.
 5. An explicit intermediate aggregation scheme for MPI-IO is implemented to mitigate the IO burden with a massive number of cores.

# Authors
- Ki-Ha Kim (k-kiha@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
- Xiaomin Pan (sanhepanxiaomin@gmail.com), Department of Mathematics, Shanghai University
- Jung-Il Choi (jic@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University

# Cite
Please use the following BibTeX when you refer to this project.

    @misc{PaScaL_TCS2022,
        title  = "Parallel and Scalable Library for Turbulent Convection Solver",
        author = "Kim, Ki-Ha and Kang, Ji-Hoon and Choi, Jung-Il",
        url    = "https://github.com/MPMC-Lab/PaScaL_TCS",
        year   = "2022"
    }


# References
For more information, please see the reference paper and [Multi-Physics Modeling and Computation Lab.](https://mpmc.yonsei.ac.kr/)

# Information
- problem: 0. RBC 1. Channel 2. Urban
- LES: 0. No-option 1. Constant Smagorinsky 2. Constant Vreman
