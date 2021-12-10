def find_zams(logl,loglh,model):
    zams=1
    while (loglh[zams] < 1.0*logl[zams]): 
     zams=zams+1
    return zams; 

def find_h(dh,center_h1,model):
    zams=1
    while (center_h1[zams] > (center_h1[1] - dh)): 
     zams=zams+1
    return zams; 

def find_mams(center_h1,model):
    mams=1
    while (center_h1[mams] > 0.5 * center_h1[1]): 
     mams=mams+1
    return mams; 

def find_tams(center_h1,model):
    tams=1
    while (center_h1[tams] > 0.05): 
     tams=tams+1
    return tams;    

def find_max(a,b,c,d):
    z= [0] * len(a)
    for i in range(0, len(a)):
      z[i]=max(a[i],b[i],c[i],d[i])   
    return z;

def binary(z):  
    for i in range(0, len(z)):
      if z[i] > 0:
         z[i] = 1
      else: 
         z[i] = 0
    return z;

def pbeta(p,b):    # Plasma Beta
    beta= 8*3.1415*p/(b**2.0)
    return beta;

def rossby(vcon,vroteq,hp,req): # Rossby number, defined as Prot/Pcon
    ross =(3.14*req/hp)*(vcon/veq)
    return ross; 
def beq(vcon,rho): # Calculate equipartition Bfield
    b=vcon*(4*3.14*rho)**0.5
    return b;

def beq(vcon,rho): # Calculate equipartition Bfield
    b=vcon*(4*3.14*rho)**0.5
    return b;

def beq_omega(vcon,rho,hp,omega): # Calculate rotating equipartition field (See eq.15 in JC20)
    b=vcon*(4*3.14*rho)**0.5
    b=b*(1+(omega*hp/vcon))**0.5
    return b;

  
# Find middle of main sequence (in age)

def find_mid_ms(model,star_age,zams,tams):
    mid_ms=1
    age_ms=(star_age[tams]-star_age[zams])/2
    while (star_age[mid_ms] < age_ms): 
     mid_ms=mid_ms+1
    return mid_ms;  

def find_frac_ms(model,star_age,zams,tams,frac):
    frac_ms=1
    age_frac_ms=(star_age[tams]-star_age[zams])*frac
    while ((star_age[frac_ms] - star_age[zams]) < age_frac_ms): 
     frac_ms=frac_ms+1
    return frac_ms;   