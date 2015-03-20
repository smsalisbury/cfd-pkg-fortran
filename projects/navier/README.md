# 2D Navier-Stokes Solver

######  Author: Spencer Salisbury
<hr />

### Description
This project solves the Navier-Stokes equations for a specified 2D enclosure. This project was developed in March and April of 2015.

Check the repository README.md for general information, such as configuration and compilation instructions.

### Usage
This project uses namelists to allow varying multiple input values without recompiling the project. The namelist variables must be in a file called "inputs". Below is a list of the namelist values and their descriptions:

<table>
<tr><th>Variable</th><th>Description</th><th>Required</th><th>Default value</th></tr>

<tr><th colspan=99>FLUID_PROPERTIES</th></tr>
<tr><td>density</td><td>The density of the fluid.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>viscosity</td><td>The viscosity of the fluid.</td><td>Yes</td><td>N/A</td></tr>

<tr><th colspan=99>BOUNDARY_CONDITIONS</th></tr>
<tr><td>bc_top_type</td><td>The boundary type of the top boundary. See "Boundary Conditions" for options.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>bc_right_type</td><td>The boundary type of the right boundary. See "Boundary Conditions" for options.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>bc_bottom_type</td><td>The boundary type of the bottom boundary. See "Boundary Conditions" for options.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>bc_left_type</td><td>The boundary type of the left boundary. See "Boundary Conditions" for options.</td><td>Yes</td><td>N/A</td></tr>

<tr><th colspan=99>INITIAL_CONDITIONS</th></tr>

<tr><th colspan=99>MESH_PROPERTIES</th></tr>

<tr><th colspan=99>ITERATIVE_PROPERTIES</th></tr>
</table>

##### Boundary Conditions
