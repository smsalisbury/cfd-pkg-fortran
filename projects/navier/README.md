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
<tr><td>bc_top_type</td><td>The boundary type of the top boundary. See "Boundary Conditions" for options.</td><td>No</td><td>wall</td></tr>
<tr><td>bc_right_type</td><td>The boundary type of the right boundary. See "Boundary Conditions" for options.</td><td>No</td><td>wall</td></tr>
<tr><td>bc_bottom_type</td><td>The boundary type of the bottom boundary. See "Boundary Conditions" for options.</td><td>No</td><td>wall</td></tr>
<tr><td>bc_left_type</td><td>The boundary type of the left boundary. See "Boundary Conditions" for options.</td><td>No</td><td>wall</td></tr>
<tr><td>bc_top_u_velocity</td><td>The u velocity at the top wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_top_v_velocity</td><td>The v velocity at the top wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_right_u_velocity</td><td>The u velocity at the right wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_right_v_velocity</td><td>The v velocity at the right wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_bottom_u_velocity</td><td>The u velocity at the bottom wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_bottom_v_velocity</td><td>The v velocity at the bottom wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_left_u_velocity</td><td>The u velocity at the left wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>
<tr><td>bc_left_v_velocity</td><td>The v velocity at the left wall. Used for 'velocity' type boundaries.</td><td>No</td><td>0</td></tr>

<tr><th colspan=99>INITIAL_CONDITIONS</th></tr>
<tr><td>ic_type</td><td>The type of initial condition to use. See "Initial Conditions" for options.</td><td>No</td><td>uniform</td></tr>

<tr><th colspan=99>MESH_PROPERTIES</th></tr>
<tr><td>x_steps</td><td>The number of steps in the x-direction to use.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>y_steps</td><td>The number of steps in the y-direction to use.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>x_size</td><td>The total x-dimension size of the mesh.</td><td>Yes</td><td>N/A</td></tr>
<tr><td>y_size</td><td>The total y-dimension size of the mesh.</td><td>Yes</td><td>N/A</td></tr>

<tr><th colspan=99>ITERATIVE_PROPERTIES</th></tr>
<tr><td>relax_mom</td><td>Momentum relaxation factor.</td><td>No</td><td>0.5</td></tr>
<tr><td>relax_press</td><td>Pressure relaxation factor.</td><td>No</td><td>0.3</td></tr>
<tr><td>relax_press_cor</td><td>Pressure correction relaxation factor.</td><td>No</td><td>1.2</td></tr>
<tr><td>conv_error</td><td>Convergence condition for continuity.</td><td>No</td><td>1.0E-5</td></tr>
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
