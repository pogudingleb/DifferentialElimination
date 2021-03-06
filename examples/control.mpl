# Example 4 from the paper https://arxiv.org/abs/1610.04022 
# The system is of the form usual in control theory, where x1, x2, x3 are state variables and u1, u2 are controls

read "../DifferentialElimination.mpl":

vars := [x1_, x2_, x3_, u1_, u2_];
eqs := [
  x1_1 - u1_0 - x2_0,
  x2_1 - u2_0 - x1_0,
  x3_1 - x3_0 * (u1_0 + u2_1 - x3_0)
]:
print(ConvertToDiffForm(eqs, vars));

for i from 1 to 3 do
  print(cat("Eliminating u1_, u2_, and ", vars[i])):
  to_eliminate := [vars[i], u1_, u2_]:
  if i = 3 then
    eqs := [op(eqs), Differentiate(eqs[2], vars, 1)]:
  end if:
  print(cat("The bound given by Theorem 2.3 is ", ComputeBound(eqs, to_eliminate))):
  prolongations := CheckPossibilityElimination(eqs, vars, to_eliminate):
  if prolongations = "no elimination" then
    print("Elimination is not possible with probability at least 99%"):
  else
    print(cat("Elimination is possible after ", prolongations, " prolongations")):
    print(cat("The result is ", ConvertToDiffForm(DifferentialElimination(eqs, vars, to_eliminate, prolongations), vars))):
  end if:
end do:

