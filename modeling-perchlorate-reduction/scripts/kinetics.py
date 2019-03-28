# Monod equation
# https://en.wikipedia.org/wiki/Monod_equation
def mm(x, km):    
    return x / (km + x)

# Monod equation altered with competitive inhibition
# as in Dudley et. al 2008
def ci(ss,ks):
    
    # ss = Substrate concentration(s), perchlorate and/or chlorate
    # ks = Substrate affinity(ies), perchlorate and/or chlorate
    
    # Allow lists to be one item
    if type(ss) is not list: 
        ss = [ss] 
    if type(ks) is not list: 
        ks = [ks] 

    if len(ss) > 1:
        # Competition for enyzme
        equation = ss[0] / (ss[0] + ks[0] * (1 + ss[1] / ks[1]))
    else:
        # No competition for enyzme
        equation = ss[0] / (ss[0] + ks[0])
            
    return equation


