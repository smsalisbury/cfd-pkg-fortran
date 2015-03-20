# 2D Navier-Stokes Solver

######  Author: Spencer Salisbury
-------------------------------

### Description
This project solves the Navier-Stokes equations for a specified 2D enclosure. This project was developed in March and April of 2015. Note that for now, this project only works for laminar flows.

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
The following are possible boundary types:
* wall
* velocity

###### Wall
The "wall" boundary type is straight-forward. It sets the velocities in both directions to 0 at the boundaries. No additional specifications are necessary.

###### Velocity
The "velocity" boundary type is used if a uniform velocity is to be set across the boundary. It is useful for problems such as the classic Lid Driven Cavity problem. If this type is chosen, the corresponding boundary velocitys should also be set. If it is not set, the values default to 0.
<table>
<tr><th>Boundary</th><th>Velocity variables</th></tr>
<tr><td rowspan=2>Top</td><td>bc_top_u_velocity</td></tr>
<tr><td>bc_top_v_velocity</td></tr>
<tr><td rowspan=2>Right</td><td>bc_bottom_u_velocity</td></tr>
<tr><td>bc_bottom_v_velocity</td></tr>
<tr><td rowspan=2>Bottom</td><td>bc_right_u_velocity</td></tr>
<tr><td>bc_right_v_velocity</td></tr>
<tr><td rowspan=2>Left</td><td>bc_left_u_velocity</td></tr>
<tr><td>bc_left_v_velocity</td></tr>
</table>

##### Initial Conditions

##### Mesh Properties

##### Iterative Properties
