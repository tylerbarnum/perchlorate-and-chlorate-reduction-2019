# Calculation of Gibbs free energy in electron (e-) equivalents
def redox_to_Ga(n,E):
    
    # E: V, redox potential
    # n: number of electrons
    
    F = 96.48 # kJ/V Faraday constant
    Ga = (- n * F * E) / n # kJ/e- Energy of electron acceptor half-reaction, normalized
    
    return Ga

# Calculation of fractional carbon/electron donor use to cell synthesis or energy production
def energy_to_fractions(Ga):

    # Ga: kJ/e- Energy of electron acceptor half-reaction

    Gc = 27.4 # kJ/e- Energy of acetate (carbon source)
    Gpyr = 35.09 # kJ/e- Energy of pyruvate (representative intermediate)
    Gp = Gpyr - Gc # kJ/e- Energy of C to pyruvate
    
    if Gp > 0:
        n = 1
    else:
        n = -1    
    
    Gpc = 18.8 # kJ/e- Energy of pyruvate and ammonium to cellular carbon
    e = 0.35 # Transfer efficiency (theoretical)
    Gd = 27.40 # kJ/e- Energy of acetate (electron donor)
    Gr = Ga - Gd # Energy per equiv. oxidized for energy production
    A = - ((Gp / e**n) + (Gpc / e)) / (e * Gr) # Equiv. of donor used for energy per cells formed
    fs = 1 / (1 + A) # Fraction of donor used for synthesis (unitless)
    fe = 1 - fs # Fraction of donor used for energy (unitless)
    
    return Gr, A, fs, fe
