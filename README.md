# bakalarkaFY

Currently, this project is attempting to solve the following two problems using the Monte-Carlo approach.

* `2reaction` (problem specified here: http://www.zdplaskin.laplace.univ-tlse.fr/index.html@p=292.html) 
* `MicroCathode` (problem specified here: http://www.zdplaskin.laplace.univ-tlse.fr/index.html@p=310.html)

Files related specifically to any of these projects are in the corresponding project folder.

# How the solver works

## Input files
First, you need to specify the parameters and the reactions in the input file (extension `.input`). Example input file:

```
time_ini = 0
time_end = 3.0e-7
calc_step = 1

# gas temperature, K
gas_temperature = 301.0

# reduced electric field, Td
reduced_field = 50.0

e = 1
Ar^+ = 1
Ar = 2.5e7

# rates in micrometers
Ar + e => e + e + Ar^+      ! 0.487e1
e + Ar^+ + Ar => Ar + Ar    ! 1.0e-1
```

Empty lines and lines starting with the hashtag sign '#' are ignored. Lines beginning with '#' can be used as comments.
The remaining lines must either set up a parameter or specify a reaction. 
There must be only one parameter or reaction specification per line
and each parameter or reaction specification must happen on one line only.

### Parameters
**Parameters** are set using the equals sign '=': `<parameter name> = <parameter value>`. Whitespace is ignored.

Several parameters must be specified, and the parser will raise an exception if
they are missing from the input file. These are:

* `time_end` - the ending time of the simulation

Most parameters are necessary for the computation, but if they are missing from the input file, 
they are set to their default value. The parser will still raise a warning to let the user know, 
that the parameter was set to its default value internally. These parameters are:

* `time_ini` - the initial time of the computation (default value 0)
* `calc_step` - specifies the number of steps of the calculation (number of loops
of the Monte-Carlo algorithm) after which the concentration values of all the species 
are saved (default value 1)

Aside from these, each species present in the reactions must have specified its initial density
(`<species name> = <initial density>`), if it is not specified, the parser will raise a warning
and set the initial density to 0.

### Reactions
**Reactions** are specified using the arrow operator '=>': `<reactants> => <products> ! <rate specification>`.

Reactants and products are separated by the plus sign '+' with a space on each side of '+'. 
The plus sign can be omitted entirely; the reactants and products can be separated by just whitespace instead.

`<rate_specification>` can either be a constant or point to a table with the rate specification.
In the example input file above, all rates are set as constants. If you wish to specify
the rates inside a table, write `table: <table name>` after the exclamation mark '!'.
Example reaction with table rate: `e + Ar => e + e + Ar^+         !   table: Ar -> Ar^+`. 
In the example, the table's name is 'Ar -> Ar^+' (surrounding whitespace is ignored).

:warning: **Note**: If you use table rate specifications, you must set the parameter
`table_file` to point to the table file containing all specified tables. 

## Table files
Table files have a much simpler syntax. 
Empty lines and lines starting with the hashtag sign '#' are ignored.

When you want to specify a table in the file, you must start by writing the table name
(the first line of the table). On the following lines, you can write any comments, and these lines are ignored.
You do not need to use the hash sign here. When you are done with commenting,
insert a dashed line ('-----'). It should contain at least 2 dashes.

After the dashed line, each line is considered a table row. Each table row must contain
exactly **two columns**, these are parsed and saved as floats. 

The table ends with another dashed line ('-----'). Again, at least 2 dashes should be used.

The parser will then linearly interpolate the data from the table.

Example table file:
```
Ar -> Ar^+
  1.580e+1 / threshold energy
COMMENT: RAPP-SCHRAM
UPDATED: 2010-03-02 16:19:07
------------------------------------------------------------
  1.580e+1   0.000e+0
  1.600e+1  2.020e-22
  1.700e+1  1.340e-21
  1.000e+4  1.350e-21
------------------------------------------------------------

# HERE COMES ANOTHER TABLE

Ar* -> Ar^+
   4.30  / threshold energy
--------------------------------
0.000  0.000
4.80   0.000
5.00   0.2020E-21
--------------------------------
```

## Running the solver

Once you have the input file (and table file if necessary) ready, you can run
the algorithm.

### Method `parse_file` from `input_parse.py`
First, you need to parse the input file using the function `parse_file` from `input_parser.py`:

`all_species, parameters, reactions = parse_file(input_filename)`,

where `input_filename` is a path to the input file. Variable `all_species` 
is a set of all species specified in the reactions. Variable `parameters` is
a dictionary storing all parameters (either set in the input file or the parser to the default values). Variable `reactions` is a list of objects of type `Reaction` storing all necessary information about each reaction specified in the input file.

Now, if necessary, you can set more parameters using regular Python 3 syntax:

`parameters[parameter_name] = parameter_value`

### Method `solve` from `solver.py`

Suppose you have a custom parameter that must be updated during the run of the
algorithm (e.g., when the reaction rates depend on it). In that case, you can create a
custom function taking `parameters` as input, modifying them, and returning
the modified dictionary. This function can then be passed into the 
`solve` method using its `update` parameter. The default value for `update` is
`None` - meaning that no custom parameters are updated during the algorithm.
If specified, the function passed in `update` is called at the end of each
loop to update the parameters.

In our example, we do not change or modify the parameters, so we call the 
`solve` method directly:

`times, values = solve(all_species, parameters, reactions)`.

The variable `times` is a list of time values at the start of each (`calc_step`-th) loop
in the algorithm. The variable `values` is a dictionary containing lists of concentration values for each specified species at times specified in the variable
`times` (in the same order).  

### Method `plot` from `plot.py`

Now, we can plot the output values of the Monte-Carlo calculation using the
`plot` method from `plot.py`:

`plot(times, values, all_species)`.

The third argument of the method specifies the species that we want to plot.
In this example, we want to plot all three species.

# Missing documentation topics:

* non-constant parameter example
* specifying custom rates dependant on custom parameters
* project file description 

