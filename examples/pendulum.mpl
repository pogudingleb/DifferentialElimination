# Example 3 from the paper https://arxiv.org/abs/1610.04022
# Pendulum from Gear and Petzold "ODE methods for the solution of differential/algebraic systems" (p. 275) with added extrnal force

read "../DifferentialElimination.mpl":

vars := [
  x, y,    # coordinates
  T,       # string tension 
  F1_, F2_ # external force
];

eqs := [
  x2 - T0 * x0 - F1_0,
  y2 + T0 * y0 - 1 - F2_0,
  x0^2 + y0^2 - 1,
  Differentiate(x0^2 + y0^2 - 1, vars, 5),
  Differentiate(x0^2 + y0^2 - 1, vars, 5, 2)
]:

print(ConvertToDiffForm(eqs, vars));

vars_to_eliminate := [
  [y, T, F1_, F2_],
  [x, T, F1_, F2_],
  [y, T, F2_],
  [x, T, F1_]
];

for i from 1 to nops(vars_to_eliminate) do
  B := ComputeBound(eqs, vars_to_eliminate[i]):
  print(cat("Bound is equal to ", B));
  elim := CheckPossibilityElimination(eqs, vars, vars_to_eliminate[i]):
  if elim = "no elimination" then
    print("Elimination is impossible with probebility at least 99%");
  else
    print(cat("Elimination is possible after ", elim, "prolongations"));
  end if:
end do:
