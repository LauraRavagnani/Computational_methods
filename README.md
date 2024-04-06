# Computational Methods in Physics

Welcome to the repository for the course "Computational Methods in Physics." This repository contains both weekly exercises and the final project for the course, aimed at providing hands-on experience with computational techniques in physics.

## Course Description

The course covers various computational methods commonly used in physics, including numerical integration, solving differential equations, Monte Carlo simulations, and more. 
      
## Contents

1. **Weekly Exercises**: Weekly exercises covering different topics and techniques taught during the course.
2. **Final Project**: A comprehensive final project that integrates multiple computational methods to solve a challenging physics problem.

## Weekly Exercises

- **Pendolo**: consider a pendulum of length $l$. Compute the ratio of its period and its period with small oscillation approximation, as a function of the maximum angle $\theta_0$. Used RUNGE-KUTTA propagator;
- **Pianeti**: simulation of the motion of celestial bodies of the Solar System (Sun, planets and Halley's comet) using VELOCITY-VERLET propagator;
- **Buca**: use SHOOTING METHOD and NUMEROV PROPAGATOR to compute eigenenergies and eigenstates of an electron in an infinite potential hole of length $L$: $$ -\frac{\hbar^2}{2m}\nabla^2\psi(x) + V(x)\psi(x) = E\psi(x) $$ studiyng the three cases:
    1. Flat hole
    2. Step potential
    3. Armonic potential
- **Schrödinger**: solve the 1-dim time-dependent Schrödinger equation using CRANK NICOLSON algorithm. The wave function annihilates at the borders. Initial wave function is: $$ \psi(x, t_0) = e^{iqx} \times e^{\frac{-(x-x_0)^2}{2\sigma^2}} $$


## Final Project

Simulation of the space probe Voyager 1 trajectory across the Solar System using VELOCITY-VERLET propagator. Gravity assist effects due to Jupiter and Saturn are highlighted.
