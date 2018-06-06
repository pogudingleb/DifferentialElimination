read "../DifferentialElimination.mpl":

# Unknown functions
vars := [x1_, x2_, x3_, u1_, u2_];
eqs := [
  x1_1 - x2_0 - u1_1, 
  x2_1 - x1_0 - u2_0, 
  x3_1 - x3_0 * (u1_0 + u2_0)
];

for pair in [[x1_, x2_], [x1_, x3_], [x2_, x3_]] do
  number_prolongations := CheckPossibilityElimination(eqs, vars, select(v -> not v in pair, vars)):
  if number_prolongations = -1 then
    print(cat("None equation only in ", pair[1], " and ", pair[2], " is a consquence of the system with probability at least 99%"));
  else
    print(cat("An equation only in ", pair[1], " and ", pair[2], " can be obtained after ", number_prolongations, " prolongations"));
    result := DifferentialElimination(eqs, vars, select(v -> not v in pair, vars), number_prolongations):
    print(cat("And the result is ", result));
  end if; 
end do;
