# Example 1 from the paper https://arxiv.org/abs/1610.04022
# Limiting case of the Van der Pol oscillator, taken from Kunkel and Mehrmann "Differential-Algebraic Equations: Analysis and Numerical Solution" (Example 1.7)

read "../DifferentialElimination.mpl":

vars := [y, z];
eqs := [
  y1 - z0,
  (1 - y0^2) * z0 - y0
]:

print(ConvertToDiffForm(eqs, vars));

result := DifferentialElimination(eqs, vars, [y], 1):
print(ConvertToDiffForm(result, vars));

