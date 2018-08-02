# Example 2 from the paper https://arxiv.org/abs/1610.04022
# The system describes coexistence of two species (prey and predators), also known as the Lotka-Volterra system

read "../DifferentialElimination.mpl":

vars := [
  x, # the population of prey
  y  # the population of predators
];

eqs := [
  x1 - alpha * x0 + beta * x0 * y0,
  y1 + gamm * y0 - delta * x0 * y0
]:

print(ConvertToDiffForm(eqs, vars));

B := ComputeBound(eqs, [y]):
print(cat("The bound is ", B));
elim_result := DifferentialElimination(eqs, vars, [y]):
print(cat("The result of elimination is ", ConvertToDiffForm(elim_result, vars)));
