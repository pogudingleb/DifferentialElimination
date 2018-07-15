# Implementation of algorithms based on a paper
# "Bounds for  elimination of unknowns in systems of differential-algebraic  equations"
# by A. Ovchinnikov, G. Pogudin, and N. Thieu Vo

  # Computes B(m, D) defined in equation (9) in the paper
  B := proc(dim, deg) 
  local result, i;
    result := 0:
    for i from 0 to dim do
      result := result + deg^(2^(i + 1) - 2):
    end do:
    result;
  end proc:

  # Computation corresponding to formula (17)
  ComputeBoundRaw := proc(m, d0, d1, h_sum, r)
  local i, result;
    result := 0:
    for i from 0 to m do
      result := result + B(i, d0 * d1^(h_sum - i - 1)):
    end do:
    result := result * d0 * d1^(min(h_sum, r) - 1):
    result;
  end proc:

  # depending on the value of "radical" parameter, computes a bound given by Theorem 2.1 or Theorem 2.2
  ComputeBound := proc(eqs, vars_to_eliminate, radical := true) 
  local diff_vars_to_eliminate, m, d0, d1, h, i, j, ord, h_sum, ideal, 
  components, degrees, hpoly, dim, deg, result;
    h := GetOrdersForSystem(eqs, vars_to_eliminate):
    diff_vars_to_eliminate := {}:
    for i from 1 to nops(h) do
      for j from 0 to h[i] - 1 do
        diff_vars_to_eliminate := { op(diff_vars_to_eliminate), MakeDerivative(vars_to_eliminate[i], j) };
      end do:
    end do:
    m := Groebner[HilbertDimension](eqs, diff_vars_to_eliminate):

    if radical then
      # bound given by Theorem 2.2
      ideal := PolynomialIdeals[PolynomialIdeal](eqs, variables = diff_vars_to_eliminate):
      if PolynomialIdeals[IsRadical](ideal) = false then
        print("The ideal is not actually radical!");
      end if:
      components := [PolynomialIdeals[PrimeDecomposition](ideal)]:
      degrees := [seq(0, i = 0..m)]:
      for i from 1 to nops(components) do
        hpoly := Groebner[HilbertPolynomial]( PolynomialIdeals[Generators](components[i]), tdeg(op(diff_vars_to_eliminate), aux_var), hilb_var ):
        dim := degree(hpoly, hilb_var):
        deg := coeff(hpoly, hilb_var, dim) * (dim)!:
        if hpoly <> 0 then
          degrees[dim + 1] := degrees[dim + 1] + deg:
        end if:
      end do:
      result := 0:
      for j from 0 to m do
        for i from 0 to j do
          if degrees[j + 1] > 0 then
            result := result + degrees[j + 1]^( 2 * (2^i - 1) ):
          end if:
        end do:
      end do:
      result;
    else
      # bound given by Theorem 2.1
      d0 := min( seq(degree(eqs[i], diff_vars_to_eliminate), i = 1..nops(eqs)) ):
      d1 := max( seq(degree(eqs[i], diff_vars_to_eliminate), i = 1..nops(eqs)) ):
      h_sum := foldl(`+`, op(h)):
      ComputeBoundRaw(m, d0, d1, h_sum, nops(eqs));
    end if:
  end proc:

  # Checks possibility of elimination based on the theory developed in Section 5
  # consistency of the resulting polynomial system can be checked using Groebner bases or triangular sets
  # (method = "groebner" or method = "triangular")
  CheckPossibilityElimination := proc(eqs, vars, vars_to_eliminate, p := 0.99, radical := true, method := "groebner") 
  local i, B, prolonged, h, h_sum, vars_to_keep, diff_vars_to_eliminate, diff_vars_to_keep, degrees, sample_size, roll,
  specialization, gb, result, ord, prolonged_spec, ring, tr;
    B := ComputeBound(eqs, vars_to_eliminate, radical);
    prolonged := eqs;
    h := GetOrdersForSystem(eqs, vars):
    h_sum := foldl(`+`, op(h)):
    vars_to_keep := select(v -> not v in vars_to_eliminate, vars):
    result := "no elimination":

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
      prolonged_spec := subs(specialization, prolonged);
      if method = "groebner" then
        gb := Groebner[Basis]( subs(specialization, prolonged), tdeg(op(diff_vars_to_eliminate)) ):
        if gb = [1] then
          result := i:
          break:
        end if:
      elif method = "triangular" then
        ring := RegularChains[PolynomialRing]([op(diff_vars_to_eliminate)]);
        tr := RegularChains[Triangularize](prolonged_spec, ring);
        if (nops(tr) = 0) or (RegularChains[Equations](tr[1], ring) = [1]) then
          result := i:
          break:
        end if:
      else
        print("No such method");
      end if:
    end do:

    result:
  end proc:

  # If number_prolongations = -1, then the algorithm tries all possible numbers form 1 to the bound
  # If number_prolongations => 0, then the algorithm performs the specified number of prolongations
  DifferentialElimination := proc(eqs, vars, vars_to_eliminate, number_prolongations := -1, radical := true) 
  local B, prolonged, h, i, ord, h_sum, vars_to_keep, result, diff_vars_to_eliminate, diff_vars_to_keep, eliminant;
    B := ComputeBound(eqs, vars_to_eliminate, radical);
    prolonged := eqs;

    h := GetOrdersForSystem(eqs, vars):
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
    result := [seq(-1, i = 0..nops(vars))];
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

  GetOrdersForSystem := proc(eqs, vars)
  local i, h, ord;
    h := [seq(-1, i = 1..nops(vars))]:
    for i from 1 to nops(eqs) do
      ord := GetOrders(eqs[i], vars):
      h := [seq(max(h[i], ord[i] + 1), i = 1..nops(vars))]:
    end do:
    h;
  end proc:

  ##########################################

  ConvertToDiffForm := proc(eqs, vars)
  local ords, subst, i, j, v;
    ords := GetOrdersForSystem(eqs, vars):
    subst := {}:
    for i from 1 to nops(vars) do
      v := vars[i]:
      for j from 0 to ords[i] do
        subst := {
          op(subst),
          MakeDerivative(v, j) = diff(v(t), [t$j]) 
        }:
      end do:
    end do:
    subs(subst, eqs);
  end proc:

#####################################################################


