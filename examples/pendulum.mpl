read "../DifferentialElimination.mpl":

vars := [x, y, u, v, T, F];
eqs := [
  x1 - u0,
  y1 - v0,
  u1 - T0 * x0,
  v1 + T0 * y0 - 1 - F0,
  x0^2 + y0^2 - 1,
  x0 * u0 + y0 * v0,
  Differentiate(x0 * u0 + y0 * v0, vars, 5)
]:

print(ConvertToDiffForm(eqs, vars));

for w in ["with", "without"] do
  for var in [x, y] do
    print(cat("Obtaining equation only in ", var, " ", w, " external force"));
    eqs_s := eqs:
    if w = "without" then
      eqs_s := subs({F0 = 0}, eqs_s):
    end if:
    if var = x then
      to_eliminate := [y, u, v, F, T]:
    else
      to_eliminate := [x, u, v, F, T]:
    end if:
    print(to_eliminate):
    print(cat("Bound given by the theorem ", ComputeBound(eqs_s, to_eliminate, true)));
    prolongations := CheckPossibilityElimination(eqs_s, vars, to_eliminate, 0.99, true);
    if prolongations = "no elimination" then
      print("Elimination is impossible with probability at least 99%"):
    else
      print(cat("Elimination is possible after ", prolongations, " prolongations with probability at least 99%")):
      result := DifferentialElimination(eqs_s, vars, to_eliminate, prolongations):
      print(cat("The result is ", ConvertToDiffForm(result, vars))):
    end if:
  end do:
end do:

