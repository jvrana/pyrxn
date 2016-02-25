# pyrxn
An easy to use chemical reaction network constructor for mass action kinetics.

<h1> Usage </h1>

pyrxn is super simple to use. You create mass action kinetics using strings and floats/ints as parameter values. Pyrxn will take care
of constructing the stoichiometery matrices and rate laws. This is 
particularly useful for quick exploration of different mass action models that you don't want to invest alot of time in. After the primary
exploration, using a heavy duty tool such as Mathematica or COPASI is appropriate.

<h2> Examples </h2>

First you initialize a chemical reaction network by calling CRN(). 
    
    from pyrxn import Reaction, CRN
    c = CRN()

To add reactions, you simply call .r(). You construct reactions using a string as an argument. 
Parameter values are supplied as ints/floats as arguments. For example:

    c.r("A + B > AB", 1) # an irreversible second order reaction composed of A and B with a rate constant of kf=1
    c.r("3AB + C <> X", 3, 1) # a reversible fourth order reaction with rate constants kf=3, kr=1. Rates and exponents are built automatically
    c.r("X > ", 0.1) # irreversible first order reaction
    c.r(" > R", 0.01) # irreversible zero order reaction

To run the trajectory on your reaction, you initialize your CRN with a dictionary, then call .run(step_size, t_final). Data is returned as
a pandas dataframe, ready for plotting.

    c.initialize({"A": 4, "B": 3}) # initializes A and B to 4 and 3 respectively. All other species are assumed to be 0.
    df = c.run(0.01, 100) # Run with step_size of 0.01 to a final time of 100

You can perform a dose response.
    dose = 10**np.linspace(-3, 2, 50)
    c.dose_response("C", dose, 0.01, 100) # dose response on C from 1e-3 to 1e2 with step_size 0.01 to a final time of 100 each run
    
Heres an example of a three ring oscillator with dampening oscillations:

    c = CRN()
    c.r("PA > PA + A", 1)
    c.r("PB > PB + B", 1)
    c.r("PC > PC + C", 1)
    c.r("PC + 2B <> PC*", 0.5, 0.1)
    c.r("PB + 3A <> PB*", 0.5, 0.1)
    c.r("PA + C <> PC*", 0.5, 0.1)
    c.r("A > ", 0.1)
    c.r("B > ", 0.1)
    c.r("C > ", 0.1)
    c.initialize({"PA": 1, "PB": 1, "PC": 1})
    c.run(0.01, 1000).plot()
