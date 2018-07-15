read "../DifferentialElimination.mpl":

vars := [y, z];
eqs := [
  y1 - z0,
  (1 - y0^2) * z0 - y0
]:

print(ConvertToDiffForm(eqs, vars));

result := DifferentialElimination(eqs, vars, [y], 1):
print(ConvertToDiffForm(result, vars));

