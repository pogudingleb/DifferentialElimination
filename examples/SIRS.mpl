read "../DifferentialElimination.mpl":

# Reference:
# S. Busenberg and K. Cooke, "Vertically Transmitted Diseases: Models and Dynamics", Section 2.14
# url: http://dx.doi.org/10.1007/978-3-642-75301-5

# Unknown functions
# S - Susceptible, In - infected, R - recovered
vars := [S, In, R];
# scalar parameters, in this example specialized to random numbers
parameters := [p, q, b1, b2, b3, d, e, c, eps, lambda, delta];
roll := rand(100);
subst := map(v -> v = roll(), parameters);
eqs := subs(subst, [
  S1 + In1 - b1 * S0 - p * b2 * In0 - b3 * R0 + d * S0 - e * R0 - q * b2 + (d + eps + c) * In0,
  (S0 + I0 + R0) * (In1 - q * b2 + (d + eps + c) * In0) - lambda * S0 * In0,
  R1 + (d + delta + e) * R0 - c * In0
]);

number_prolongations := CheckPossibilityElimination(eqs, vars, [In, R]):
cat("In and R can be elimimated after ", number_prolongations, " prolongations");
result := DifferentialElimination(eqs, vars, [In, R], number_prolongations):
cat("The results of elimination is ", result);
