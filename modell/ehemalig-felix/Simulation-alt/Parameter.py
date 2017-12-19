# Simulationsparameter###########################################################
# Maximale Simulationszeit in ms
TMAX = 100

# Schrittzeit in Sekundenbruchteilen
STEPTIME = 1000
# Spannungsamplitude (an der Zelle)
AMPLITUDE = 10 

freq = 20.0  # in Hz
# Zellparameter##################################################################

# Zellradius in m
RADIUS = 1.5e-6

# Widerstand der ATP-Synthase in V/ATP/s
# Annahme: Bei 100mv Protonmotive Force werden 100 ATP/s umgesetzt
RATPSYNTHASE = 1

# Anzahl ATP-Synthasen in der Membran
#lt REchnung s. Wiki
NATPSYNTHASE = 2200
#Hoechstwert Anhand der Flaeche: 400000

# CsR Parameter
PCSRCHANNEL = 1
NCSRCHANNEL = 2200
# Umweltparameter################################################################
import numpy as np
# Ausgangs-Ph-Werte
cprotIn = np.zeros(TMAX)
cprotIn[0] = 7e-14
cprotOut = 7e-14

# Ausgangsspannungen
uout = np.zeros(TMAX)
uin = np.zeros(TMAX)
uout[0] = 1.0
uin[0] = 0.8

# Naturkonstanten################################################################
RGASCONSTANT = 8.3144598
TEMPABS = 300
FFARADAY = 96485.33289

# Errechnete Zellparameter#######################################################

# Aus dem Radius errechnetes Volumen der Zelle
VOLUME = 3.14 * 4.0 * RADIUS * RADIUS * RADIUS / 3.0

# Wir rechnen in ms Schritten und damit brauchen wir einen anderen Widerstand
RATPSYNTHASE *= STEPTIME

# Kapazitaet der Zelle
MU0 = 8.854e-12
MUR = 50
DLIPID = 2e-9
A = 4 * 3.14 * RADIUS * RADIUS
CAPACITY = MU0 * MUR * A / DLIPID

t = 0

# Anzahl der Protonen in der Zelle als Array ueber die Zeit
nprotin = np.zeros(TMAX)
nprotin[t] = cprotIn[0] * VOLUME

tseries = np.arange(TMAX)

# Wir rechnen in kleinerern Schritten
freq /= STEPTIME  