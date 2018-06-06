########### Auxiliary procedures for working with differential polynomials ################
# General convention: variables are strings not ending with digits
# derivatives of a variable v are of the form v0, v1, v2, ...

  MakeDerivative := proc(var_name, der_order):
    cat(var_name, der_order):
  end proc:

  ##########################################

  DifferentiateOnce := proc(diff_poly, var_list, max_ord) local result, i, j:
    result := 0:
    for i from 1 to nops(var_list) do
      for j from 0 to max_ord do
        result := result + diff(diff_poly, MakeDerivative(var_list[i], j)) * MakeDerivative(var_list[i], j + 1):
      end do:
    end do:
    simplify(result):
  end proc:

  ##########################################

  Differentiate := proc(diff_poly, var_list, max_ords, ord := 1) local result, i;
    result := diff_poly:
    for i from 1 to ord do
      result := DifferentiateOnce(result, var_list, max_ords):
    end do:
    result:
  end proc:

  ##########################################

  # decomposes a derivative into [variable, order], if possible
  DecomposeDerivative := proc(diff_var, vars) local matched, var, v, ord;
    matched := false:
    var := 0:
    for v in vars do
      if StringTools[IsPrefix](v, diff_var) then
        matched := true:
        var := v:
      end if:
    end do:
    if matched = true then
      ord := parse(StringTools[Drop](diff_var, length(var))):
      [var, ord];
    else
      [];
    end if:
  end proc:

  #########################################

  GetOrders := proc(diff_poly, vars) 
  local result, diff_vars, i, diff_v, decomposition;
    result := [seq(0, i = 0..nops(vars))];
    diff_vars := indets(diff_poly):
    for i from 1 to nops(vars) do
      for diff_v in diff_vars do
        decomposition := DecomposeDerivative(diff_v, [vars[i]]):
        if nops(decomposition) = 2 then
          result[i] := max(result[i], decomposition[2]):
        end if:
      end do:
    end do:
    return result;
  end proc:

#####################################################################

  # the value of R_u(dim, deg)
  Ru := proc(dim, deg) 
  local result, i;
    result := 0:
    for i from 0 to dim do
      result := result + deg^(2^(i + 1) - 2):
    end do:
    result;
  end proc:

  ComputeBoundRaw := proc(m, d0, d1, h_sum, r, radical := false)
  local i, result;
    result := 0:
    for i from 0 to m do
      result := result + Ru(i, d0 * d1^(h_sum - i - 1)):
    end do:
    if radical = false then
      result := result * d0 * d1^(min(h_sum, r) - 1):
    end if:
    result;
  end proc:

  ComputeBound := proc(eqs, vars_to_eliminate, radical := false) 
  local diff_vars_to_eliminate, m, d0, d1, h, i, ord, h_sum;
    diff_vars_to_eliminate := select(
      v -> nops(DecomposeDerivative(v, vars_to_eliminate)) = 2,
      indets(eqs)
    ):
    m := Groebner[HilbertDimension](eqs, diff_vars_to_eliminate):
    d0 := min( seq(degree(eqs[i], diff_vars_to_eliminate), i = 1..nops(eqs)) ):
    d1 := max( seq(degree(eqs[i], diff_vars_to_eliminate), i = 1..nops(eqs)) ):
    h := [seq(0, i = 1..nops(vars_to_eliminate))]:
    for i from 1 to nops(eqs) do
      ord := GetOrders(eqs[i], vars_to_eliminate):
      h := [seq(max(h[i], ord[i]), i = 1..nops(vars_to_eliminate))]:
    end do:
    h_sum := foldl(`+`, op(h)):
    ComputeBoundRaw(m, d0, d1, h_sum, nops(eqs), radical);
  end proc:

  CheckPossibilityElimination := proc(eqs, vars, vars_to_eliminate, p := 0.99, radical := false) 
  local i, B, prolonged, h, h_sum, vars_to_keep, diff_vars_to_eliminate, diff_vars_to_keep, degrees, sample_size, roll,
  specialization, gb, result;
    B := ComputeBound(eqs, vars_to_eliminate, radical);
    prolonged := eqs;
    h := [seq(0, i = 1..nops(vars))]:
    for i from 1 to nops(eqs) do
      ord := GetOrders(eqs[i], vars):
      h := [seq(max(h[i], ord[i]), i = 1..nops(vars))]:
    end do:
    h_sum := foldl(`+`, op(h)):
    vars_to_keep := select(v -> not v in vars_to_eliminate, vars):
    result := -1:

    for i from 1 to B do
      prolonged := [
        op(prolonged),
        seq( Differentiate(prolonged[-j], vars, i + h_sum), j = 1..nops(eqs) )
      ]:
      diff_vars_to_eliminate := select(
        v -> nops(DecomposeDerivative(v, vars_to_eliminate)) = 2,
        indets(prolonged)
      ):
      diff_vars_to_keep := select(
        v -> nops(DecomposeDerivative(v, vars_to_keep)) = 2,
        indets(prolonged)
      ):
      degrees := map(p -> degree(p, [op(diff_vars_to_eliminate), op(diff_vars_to_keep)]), prolonged):
      degrees := sort(degrees, `>`)[1..min(nops(degrees), nops(diff_vars_to_eliminate) + nops(diff_vars_to_keep))]:

      sample_size := ceil(foldl(`*`, op(degrees)) / (1 - p)):
      roll := rand(sample_size):
      specialization := map(v -> v = roll(), diff_vars_to_keep):
      gb := Groebner[Basis]( subs(specialization, prolonged), tdeg(op(diff_vars_to_eliminate)) ):
      if gb = [1] then
        result := i:
        break:
      end if:
    end do:

    result:
  end proc:

  DifferentialElimination := proc(eqs, vars, vars_to_eliminate, number_prolongations := -1, radical := false) 
  local B;
    B := ComputeBound(eqs, vars_to_eliminate, radical);
    prolonged := eqs;

    h := [seq(0, i = 1..nops(vars))]:
    for i from 1 to nops(eqs) do
      ord := GetOrders(eqs[i], vars):
      h := [seq(max(h[i], ord[i]), i = 1..nops(vars))]:
    end do:
    h_sum := foldl(`+`, op(h)):
    vars_to_keep := select(v -> not v in vars_to_eliminate, vars):
    result := []:

    for i from 1 to B do
      prolonged := [
        op(prolonged),
        seq( Differentiate(prolonged[-j], vars, i + h_sum), j = 1..nops(eqs) )
      ]:
      if number_prolongations <> -1 and number_prolongations <> i then
        next;
      end if:
      diff_vars_to_eliminate := select(
        v -> nops(DecomposeDerivative(v, vars_to_eliminate)) = 2,
        indets(prolonged)
      ):
      diff_vars_to_keep := select(
        v -> nops(DecomposeDerivative(v, vars_to_keep)) = 2,
        indets(prolonged)
      ):

      eliminant := PolynomialIdeals[EliminationIdeal](
        PolynomialIdeals[PolynomialIdeal](prolonged, variables = {op(diff_vars_to_eliminate), op(diff_vars_to_keep)}),
        diff_vars_to_keep
      ):
      if eliminant <> PolynomialIdeals[PolynomialIdeal]([0]) then
        result := PolynomialIdeals[Generators](eliminant):
        break:
      end if:
    end do:

    result:
  end proc:
