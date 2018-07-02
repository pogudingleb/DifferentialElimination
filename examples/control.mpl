read "../DifferentialElimination.mpl":

# Unknown functions
vars := [x1_, x2_, x3_, u1_, u2_];
eqs := [
  x1_1 - u1_0 - x2_0,
  x2_1 - u2_0 - x1_0,
  x3_1 - x3_0 * (u1_0 + u2_1 - x3_0)
];

for i from 1 to 3 do
  print(cat("Eliminating u1_, u2_, and ", vars[i])):
  to_eliminate := [vars[i], u1_, u2_]:
  if i = 3 then
    eqs := [op(eqs), Differentiate(eqs[2], vars, 1)]:
  end if:
  print(cat("The bound given by Theorem 2.3 is ", ComputeBound(eqs, to_eliminate))):
  prolongations := CheckPossibilityElimination(eqs, vars, to_eliminate):
  if prolongations = "no elimination" then
    print("Elimination is not possible"):
  else
    print(cat("Elimination is possible after ", prolongations, " prolongations")):
    print(cat("The result is ", DifferentialElimination(eqs, vars, to_eliminate, prolongations))):
  end if:
end do:

