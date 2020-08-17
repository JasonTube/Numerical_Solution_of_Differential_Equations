# Problem1
Program the central FD method for the self-adjoint BVP using a uniform grid and the central FD scheme. Test your code for some cases.
Plot the computed solution and the exact solution, and the errors for a particular grid n = 80. Do the grid refinement analysis, to determine the order of accuracy of the global solution. Also try to answer the following questions:

• Can your code handle the case when a = 0 or b = 0?

• If the central finite difference scheme is used for the special equivalent differential equation and what are the advantages or disadvantages?

# Problem2
Modify the Matlab code non_tp.m to apply the finite difference method and Newton’s non-linear solver to find a numerical solution to the non-linear pendulum model.

# Problem3
Implement and compare the Gauss-Seidel, and the SOR (trying to find the best ω by testing), methods for the elliptic problem.

Do the grid refinement analysis in the infinity norm. Take the tolerance as 10−5. Does the method behave like a second order method? Compare also the number of iterations and test the optimal relaxation factor ω. Plot the solution and the error for n = 32.

Having made sure that you code is working correctly, try your code with a point source f(x, y) = δ(c − 0.5)δ(y − 0.5) and ux = −1 at x = 1, with p(x, y) = 1 and r(x, y) = 0. Note that the u(x, y) can be interpreted as the steady state temperature distribution of a room with insulated wall on three sides, a constant heat flow in from one side, and a point source (a heater for example) in the room. Note that the heat source can be expressed as f(n/2, n/2) = 1/h2 and f(i, j) = 0 for other grid points. Use the mesh and contour plots to visualize the solution for n = 36(mesh(x, y, u), contour(x, y, u, 30)).
