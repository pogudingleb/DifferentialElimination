# DifferentialElimination
Experimental implementations of algorithms based on the paper https://arxiv.org/abs/1610.04022

## Usage

Main functions implemented in "DifferentialElimination.mpl" are

### ComputeBound(eqs, vars_to_eliminate, radical)

Computes the bound for elimination of vars_to_eliminate in the system eqs of differential-algebraic equations

* **eqs** is a list of the equations

* **vars_to_eliminate** is a list of the unknowns to be eliminated

* **radical** is an optional parameter (the default value is true). If the value is true, the algorithm uses the bound for radical ideals (from Theorem 2), otherwise the algorithm uses the bound from Theorem 1. The former bound is always smaller. If radical is set to be true, but the ideal in not actually radical, the algorithm will raise an error.

### CheckPossibilityElimination(eqs, vars, vars_to_eliminate, p, radical, method)

Checks possibility of elimination of vars_to_eliminate in the system eqs of differential-algebraic equations.
If the elimination is possible, return the minimal number of sufficient prolongations, otherwise returns "no elimination".

* **eqs** is a list of the equations

* **vars** is a list of all the unknowns

* **vars_to_eliminate** is a list of the unknowns to be eliminated

* **p** is the user-specified probability of the correcteness of the result

* **radical** is an optional parameter (the default value is true). If the value is true, the algorithm uses the bound for radical ideals (from Theorem 2), otherwise the algorithm uses the bound from Theorem 1. The former bound is always smaller. If radical is set to be true, but the ideal in not actually radical, the algorithm will raise an error.

* **method** is an optional parameter (the default value is "groebner"). The parameter specifies the method used to determine the consistency of the resulting polynomial system (for the details, see Section 5 of the paper). On different examples, one of the methods can be faster than the other.

  * **groebner** - Groebner bases
  * **triangular** - Triangular decomposition (RegularChains library)
  
### DifferentialElimination(eqs, vars, vars_to_eliminate, number_prolongations, radical)

Computes the result of elimination of vars_to_eliminate in the system eqs of differential-algebraic equations.

* **eqs** is a list of the equations

* **vars** is a list of all the unknowns

* **vars_to_eliminate** is a list of the unknowns to be eliminated

* **number_prolongations** is an optional parameter (the default value is -1). 
  * **-1** the algorithm performs the number of prolongtions given by *ComputeBound* function
  * **nonegative integer** the algorithm performs the specified number of prolongations. For efficiency, this number can be precomputed using *CheckPossibilityElimination* function
  
* **radical** is an optional parameter (the default value is true), is used only if number_prolongations = -1. If the value is true, the algorithm uses the bound for radical ideals (from Theorem 2), otherwise the algorithm uses the bound from Theorem 1. The former bound is always smaller. If radical is set to be true, but the ideal in not actually radical, the algorithm will raise an error.

## Structure

* **DifferentialElimination.mpl** contains the algorithms for elimination

* **examples/** folder contains the code corresponding to the Examples 1-4 from the paper

## Support

The software is partially supported by the National Science Foundation.
