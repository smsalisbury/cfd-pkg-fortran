# FORTRAN CFD Package

######  Author: Spencer Salisbury
<hr />

### Description
This package was written for MAE 5440 (Computational Fluid Dynamics) at Utah State University.

### Compilation
These projects use makefiles to facilitate compilation. To complie a project, simply naviate to the project directory (ex: "projects/navier/") and type the command <code>make</code> in the command line:
<pre><code>$ make</code></pre>
Note that this works in both Windows and Unix/Linux environments.
To recompile, the compiled objects need to be cleaned. To remove all object files:
<pre><code>$ make clean</code></pre>
To remove all object files and executables:
<pre><code>$ make veryclean</code></pre>

### Usage
Each project has a "config.f90" file. Use this file to set the working precision. Find the line
<pre><code>integer,parameter :: wp = sp</code></pre>
and change the equality to the desired precision prior to compilation.
<table>
<tr><th>Abbr</th><th>Precision Level</th></tr>
<tr><td>sp</td><td>Single</td></tr>
<tr><td>dp</td><td>Double</td></tr>
<tr><td>ep</td><td>Extended</td></tr>
<tr><td>qp</td><td>Quadruple</td></tr>
</table>
Other usage for each project varies. Check the README.md in each project directory to get specific usage.
