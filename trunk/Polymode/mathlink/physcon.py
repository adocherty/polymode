#
# Downloaded from: http://www.hjcb.nl/python/physcon.html
# Herman J.C. Berendsen
#
# physcon Physical constants
# Note: type physcon.help() (after import physcon)

from math import pi

# define light velocity, because we need it in calculations
cloc =299792458.

# dictionary of physical constants in SI units. Values june 2005 http://www.physics.nist.gov/cuu/Constants/
# each item: [description (string), symbol (string), value (float), sd (float), relat. sd (float),
#             value(sd) unit (string), source (string)]
constant_data={'lightvel':['velocity of light in vacuum','c',cloc,0., 0.,'299 792 458(ex) m/s', 'CODATA 2002'],
     'planck':["Planck's constant",'h',6.6260693e-34,1.1e-40,1.7e-7,'6.626 0693(11) e-34 J s', 'CODATA 2002'],
     'dirac':["Dirac's constant = h/(2 pi)",'hbar',1.05457168e-34,1.8e-41,1.7e-7,'1.054 571 68(18) e-34 J s', 'CODATA 2002'],
     'permeab':['magnetic permeability of vacuum','mu_0',4e-7*pi,0.,0.,'1.256 637 061... e-6 N A^-2',''],
     'permit':['dielectric permittivity of vacuum','eps_0',1.e7/(4*pi*cloc*cloc),0.,0.,'8.854 187 817... e-12 F/m',''],
     'gravit':['Newton constant of gravitation','G',6.6742e-11, 1.0e-14,1.5e-4,'6.6742(10) e-11 m^3 kg^-1 s^-2',
               'CODATA 2002'],
     'charge-e':['elementary charge','e',1.60217653e-19,1.4e-26,8.5e-8,'1.602 176 53(14) e-19 C','CODATA 2002'],
     'mass-e':['electron mass','m_e',9.1093826e-31,1.6e-37,1.7e-7,'9.109 3826(16) e-31 kg','CODATA 2002'],
     'mass-e/u':['electron mass in u','m_e_u',5.4857990945e-4,2.4e-13,4.4e-10,'5.485 799 0945(24) u','CODATA 2002'],
     'mass-p':['proton mass','m_p',1.67262171e-27,2.9e-34,1.7e-7,'1.672 621 71(29) e-27 kg','CODATA 2002'],
     'mass-p/u':['proton mass in u','m_p_u',1.00727646688,1.3e-10,1.3e-10,'1.007 276 466 88(13) u','CODATA 2002'],
     'mass-n':['neutron mass','m_n',1.67492728e-27,2.9e-34,1.7e-7,'1.674 927 28(29) e-27 kg','CODATA 2002'],
     'mass-n/u':['neutron mass in u','m_n_u',1.00866491560,5.5e-10,5.5e-10,'1.008 664 915 60(55) u','CODATA 2002'],
     'mass-d':['deuteron mass','m_d',3.34358335e-27,5.7e-34,1.7e-7,'3.343 583 35(57) e-27 kg','CODATA 2002'],
     'mass-d/u':['deuteron mass in u','m_d_u',2.01355321270,3.5e-10,1.7e-10,'2.013 553 212 70(35) u','CODATA 2002'],
     'mass-mu':['muon mass','m_m',1.88353140e-28,3.3e-35,1.7e-7,'1.883 531 40(33) e-28 kg','CODATA 2002'],
     'mass-mu/u':['muon mass in u','m_m_u',0.1134289264,3.0e-9,2.6e-8,'0.113 428 9264(30) u','CODATA 2002'],
     'ratio-me/mp':['electron/proton mass ratio','ratio_memp',5.4461702173e-4,2.5e-13,4.6e-10,'5.446 170 2173(25) e-4','CODATA 2002'],
     'ratio-mp/me':['proton/electron mass ratio','ratio_mpme',1836.15267261,8.5e-7,4.6e-10,'1836.152 672 61(85)','CODATA 2002'],
     'amu':['unified atomic mass unit = 1/12 m(12C)','u',1.66053886e-27,2.8e-34,1.7e-7,'1.660 538 86(28) e-27 kg','CODATA 2002'],
     'avogadro':['Avogadro constant','N_A',6.0221415e23,1.0e17,1.7e-7,'6.022 1415(10) e23 mol^-1','CODATA 2002'],
     'boltzmann':['Boltzmann constant','k_B',1.3806505e-23,2.4e-29,1.8e-6,'1.380 6505(24) J/K','CODATA 2002'],
     'gas':['molar gas constant = N_A k_B','R',8.314472,1.5e-5,1.7e-6,'8.314 472(15) J mol^-1 K^-1','CODATA 2002'],
     'faraday':['Faraday constant - N_A e','F',96485.3383,8.3e-3,8.6e-8,'96 485.3383(83) C/mol','CODATA 2002'],
     'bohrradius':['Bohr radius = 4 pi eps_0 hbar^2/(m_e e^2)','a_0',5.291772108e-11,1.8e-19,3.3e-9,'0.529 177 2108(18) e-10 m','CODATA 2002'],
     'magflux-qu':['magnetic flux quantum = h/(2 e)','Phi_0',2.06783372e-15,1.8e-22,8.5e-8,'2.067 833 72(18) Wb','CODATA 2002'],
     'conduct-qu':['conductance quantum = 2 e^2/h','G_0',7.7748091733e-5,2.6e-13,3.3e-9,'7.748 091 733(26) e-5 S','CODATA 2002'],
     'josephson':['Josephson constant = 2 e/h','K_J',4.83597879e14, 4.1e7,8.5e-8,'4.835 978 79(41) e14 Hz/V','CODATA 2002'],
     'bohrmagn':['Bohr magneton = e hbar/(2 m_e)','mu_B',9.27400949e-24,8.0e-31,8.6e-8,'9.274 009 49(80) e-24 J/T','CODATA 2002'],
     'nuclmagn':['nuclear magneton = e hbar/(2 m_p)','mu_N',5.05078343e-27,4.3e-34,8.6e-8,'5.050 783 43(43) e-27 J/T','CODATA 2002'],
     'magnmom-e':['electron magnetic moment','mu_e',-9.28476412e-24,8.0e-31,8.6e-8,'-9.284 764 12(80) e-24 J/T','CODATA 2002'],
     'magnmom-p':['proton magnetic moment','mu_p',1.41060671e-26,1.2e-33,8.7e-8,'1.410 606 71(12) e-26 J/T','CODATA 2002'],
     'gfactor-e':['electron g-factor','g_e',-2.0023193043718,7.5e-12,3.8e-12,'-2.002 319 304 3718(75)','CODATA 2002'],
     'gfactor-p':['proton g-factor','g_p',5.585694701, 5.6e-8,1.0e-8,'5.585 694 701(56)','CODATA 2002'],
     'alpha':['fine-structure constant = e^2/(4 pi eps_0 hbar c)','alpha',7.297352568e-3,2.4e-11,3.3e-9,'7.297 352 568(24) e-3','CODATA 2002'],
     'alpha-1':['inverse fine-structure constant = 4 pi eps_0 hbar c/e^2','',137.03599911,4.6e-7,3.3e-9,'137.035 999 11(46)','CODATA 2002'],
     'gyromagratio-p':['proton gyromagnetic ratio','gamma_p',2.67522205e8,23.,8.6e-8,'2.675 222 05(23) e8 s^-1 T^-1','CODATA 2002'],
     'magres-p':['magnetic resonance frequency proton = gamma_p/(2*pi)','',4.25774813e7,3.7,8.6e-8,'42.577 4813(37) MHz/T','CODATA 2002'],
     'rydberg':['Rydberg constant = alpha^2 m_e c/(2 h)','R_infty',10973731.568525,7.3e-5,6.6e-12,'10 973 731.568 525(73) m^-1','CODATA 2002'],
     'stefan-boltzm':['Stefan-Boltzmann constant = pi^2 k^4/(60 hbar^3 c^2)','sigma',5.670400e-8,4.0e-13,7.0e-6,'5.670 400(40) e-8 W m^-2 K^-4','CODATA 2002']}             
     

# many common values are also available as global constants:
global alpha,a_0,c,e,eps_0,F,G,g_e,g_p,gamma_p,h,hbar,k_B,eV
global m_d,m_e,m_n,m_p,mu_B,mu_e,mu_N,mu_p,mu_0,N_A,R,sigma,u
alpha = constant_data['alpha'][2]
a_0 =  constant_data['bohrradius'][2]
c = cloc
e =  constant_data['charge-e'][2]
eps_0 =  constant_data['permit'][2]
F =  constant_data['faraday'][2]
G =  constant_data['gravit'][2]
g_e =  constant_data['gfactor-e'][2]
g_p =  constant_data['gfactor-p'][2]
gamma_p =  constant_data['gyromagratio-p'][2]
h = constant_data['planck'][2]
hbar = constant_data['dirac'][2]
k_B =  constant_data['boltzmann'][2]
m_d =  constant_data['mass-d'][2]
m_e =  constant_data['mass-e'][2]
m_n =  constant_data['mass-n'][2]
m_p =  constant_data['mass-n'][2]
mu_B =  constant_data['bohrmagn'][2]
mu_e =  constant_data['magnmom-e'][2]
mu_N =  constant_data['nuclmagn'][2]
mu_p =  constant_data['magnmom-p'][2]
mu_0 =  constant_data['permeab'][2]
N_A =  constant_data['avogadro'][2]
R =  constant_data['gas'][2]
sigma =  constant_data['stefan-boltzm'][2]
u =  constant_data['amu'][2]
ev = constant_data['charge-e'][2]

def help():
    print 'Available functions:'
    print '[note: key must be a string, within quotes!]' 
    print '  value(key) returns value (float)'
    print '  sd(key)    returns standard deviation (float)'
    print '  relsd(key) returns relative standard deviation (float)'
    print '  descr(key) prints description with units\n'
    print 'Available global variables:'
    print '  alpha, a_0, c, e, eps_0, F, G, g_e, g_p, gamma_p, h, hbar, k_B'
    print '  m_d, m_e, m_n, m_p, mu_B, mu_e, mu_N, mu_p, mu_0, N_A, R, sigma, u\n'
    allkeys=constant_data.keys()
    allkeys.sort()
    print 'Available keys:'
    print allkeys

def value(key):
    return constant_data[key][2]

def sd(key):
    return constant_data[key][3]

def relsd(key):
    return constant_data[key][4]

def descr(key):
    print 'Description of ',key,':'
    print '  Name:               ',constant_data[key][0]
    print '  Symbol (if avail.): ',constant_data[key][1]
    print '  Value:              ',constant_data[key][2]
    print '  Standard deviation: ',constant_data[key][3] 
    print '  Relative stdev:     ',constant_data[key][4]
    print '  value(sd) unit:     ',constant_data[key][5]
    print '  Source:             ',constant_data[key][6],'\n'

