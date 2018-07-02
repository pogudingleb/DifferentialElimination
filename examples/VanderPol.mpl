read "../DifferentialElimination.mpl":

vars := [y, z];
# scalar parameters, in this example specialized to random numbers
eqs := [
  y1 - z0,
  (1 - y0^2) * z0 - y0
];

result := DifferentialElimination(eqs, vars, [y], 1):
print(result);

