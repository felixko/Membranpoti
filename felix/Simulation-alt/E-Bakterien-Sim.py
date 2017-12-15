# TODO
# Bei der veränderung der Spannung innen ist glaube ich die Elektronenladung
# nicht berücksichtigt



from math import *
from matplotlib import pyplot as plt
import numpy as np


# Simulationsparameter##########################################################

# Maximale Simulationszeit in ms
TMAX = 1000000

# Schrittzeit in Sekundenbruchteilen
STEPTIME = 1
# Spannungsamplitude (an der Zelle)
AMPLITUDE = 1
OFFSET = 1

#Maximaler PH Unterschied
PHDMAX = 0.2

freq = 0.001  # in Hz

# Zellparameter#################################################################

# Zellradius in m
RADIUS = 1.5e-6

# Widerstand der ATP-Synthase in V/ATP/s
# Annahme: Bei 100mv Protonmotive Force werden 100 ATP/s umgesetzt
# Also 300 Protonen/s
# Bei einem Volt 3000 Protonen/s
RATPSYNTHASE = 1.0 / 3000.0

# Anzahl ATP-Synthasen in der Membran
#lt REchnung s. Wiki
# 2 nullen mehr damit es sich schneller anpassen
# Ist eh die Frage was die Praxiswerte sind
# Urspruenglich 2200
NATPSYNTHASE = 22000000
#Hoechstwert Anhand der Flaeche: 400000

# CsR Parameter
PCSRCHANNEL = 0.2
RCSRCHANNEL = 0.8 / 3000.0
# Hier genauso 2 Nullen mehr
NCSRCHANNEL = 220000



# Umweltparameter###############################################################


PHSTD = 7

# Ausgangs-Ph-Werte
cprotIn = np.zeros(TMAX)
cprotIn[0] = 7e-14



cprotOut = np.zeros(TMAX)
cprotOut[0] = PHSTD * 1e-14

# Ausgangsspannungen
uout = np.zeros(TMAX)
uin = np.zeros(TMAX)
uout.fill(1)
uin[0] = 0.8



# Naturkonstanten###############################################################
RGASCONSTANT = 8.3144598
TEMPABS = 300
FFARADAY = 96485.33289
ELEMENTARYCHARGE = 1.6e-19
# Errechnete Zellparameter######################################################

# Aus dem Radius errechnetes Volumen der Zelle
VOLUME = 3.14 * 4.0 * RADIUS * RADIUS * RADIUS / 3.0

RATPSYNTHASE /= NATPSYNTHASE
# Wir rechnen in ms Schritten und damit brauchen wir einen anderen Widerstand
RATPSYNTHASE *= STEPTIME

#Das gleiche mit den Werten des Channels
RCSRCHANNEL /= NCSRCHANNEL
RCSRCHANNEL *= STEPTIME

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


# Als Aufzeichnungsvariable
pmftseries = np.zeros(TMAX)

################################################################################
# Berechnet die Protonmotiveforce anhand des Protonen- und Spannungsunterschieds
def calcProtonmotiveforce():
	# Teilt sich auf in einen Chemischen und Elektrischen Anteil

	# Wir veraendern primaer den Elektrischen Anteil, aber muessen
	# natuerlich auch den Chemischen mit beruecksichtigen

	# Fuer beide Ansaetze siehe Wikipedia: protonmotive Force

	# 1.Elektrischer Anteil:
	# Formeln:
	# dQ = r*dU
	# dU' = C*dQ*dU

	# Ich glaube hier stimmt etwas nicht
	dU = uout[t] - uin[t]
	# je geringer die Protonenkonzentration ist, desto hoeher ist der pH wert
	# d.h im inneren ist ein hoeherer pH wert als aussen
	# dh. wir ziehen das aeussere vom inneren ab, um eine positive Zahl zu haben
	dpH = 0
	global cprotIn
	global cprotOut
	# maximale Werte pruefen
	

	# 2. Chemischer Anteil
	dpH = log10(cprotIn[t]) - log10(cprotOut[t])
	# ph ist -log(cH3O) also andersherum

	pmf = dU + 2.3 * RGASCONSTANT * TEMPABS * dpH / FFARADAY
	# Ist die pmf positiv, so werden die Protonen in die Zelle gedrueckt

	if pmf != 0:
		return pmf
	else:
		return 1e-14


def boundrychecker():
	global nprotin
	global cprotOut
	#Konzentrationen anhand PH Wert
	if cprotIn[t] < 1e-14:
		cprotIn[t] = 1e-14
	elif cprotIn[t] > 14e-14:
		cprotIn[t] = 14e-14
	if cprotOut[t] < 1e-14:
		cprotOut[t] = 1e-14
	elif cprotOut[t] > 14e-14:
		cprotOut[t] = 14e-14

	#Anzahl muss positiv sein
	if nprotin[t] <= 0:
		nprotin[t] = VOLUME * 1e-14
		cprotIn[t] = 1e-14


# Berechnet die durch die ATP-Synthase bewegten Protonen
def atpMovement():
	# Durch die ATP-Synthase werden Protonen nach innen/aussen bewegt
	boundrychecker()
	iprotmoved = pmftseries[t] / RATPSYNTHASE

	# Integral der Elektronenflussstärke
	nprotmoved = iprotmoved/STEPTIME

	return nprotmoved


# Berechnet die durch den CrS Channel bewegten Protonen
def CsRMovement():
	boundrychecker()

	# Das Ganze stimmt nicht so recht
	# Wir müssen eine Funktion finden die die bewegten Protonen gut simuliert
	
	# Idee: durch CsR wird eine Kraft zusätzlich zur pmf erzeugt
	# CsR hat einen Widerstand, die fließende Ladung wird berechnet aus
	# den beiden Kräften durch den Widerstand

	# #Die erzeugte Spannung ist unabhängig von der Anzahl, der Widerstand schon

	fCsR = pmftseries[t] - PCSRCHANNEL
	nprotCsR = fCsR / RCSRCHANNEL

	if nprotCsR > 0:
		nprotCsR = 0

	# Integral über die Zeit
	# Wir teilen dadurch da die Zeit 1/Steptime ist
	nprotCsR /= STEPTIME
	


	#nprotCsR = calcProtonmotiveforce() / RATPSYNTHASE
	return nprotCsR
	#return 0


# Berechnet die Anliegende Spannung
def calcUout(tloc):
	# Es wird eine einfache sinusspannung angelegt
	uout[tloc] = sin(t*2*pi*freq)*AMPLITUDE + OFFSET

# Berechnet die Aenderung des pH-Werts in der Lösung
def calcpHout(tloc):
	cprotOut[tloc] = sin(t*2*pi*freq)*PHDMAX + PHSTD
	cprotOut[tloc] = cprotOut[tloc] * 1e-14




################################################################################
while t < TMAX-1:
	boundrychecker()
	# pmf fuer den Zeitschritt berechnen
	# Spart Rechenzeit
	pmftseries[t] = calcProtonmotiveforce()

	# Spannung aussen berechnen
	# calcUout(t)
	# calcpHout(t)
	# Alles fuer den naechsten Zeitschritt berechnen
	# Durch ATP-Synthase bewegte Protonen
	nprotin[t+1] = nprotin[t] + atpMovement()
	# Durch csR bewegte Protonen
	#nprotin[t+1] += CsRMovement()
	
	# Aus der Anzahl gewanderter Protonen die Spannungsaenderung errechnen
	# dProt = nprotin[t]-nprotin[t+1]
	# Wie veraendert die Protonenwanderung die Spannung?
	uin[t+1] = uin[t] + ((nprotin[t+1] - nprotin[t]) * CAPACITY)
	
	cprotIn[t] = nprotin[t] / VOLUME
	
	t += 1

#Letzten Zeitschritt berechnen
nprotin[t] -= atpMovement()
cprotIn[t] = nprotin[t] /VOLUME
pmftseries[t] = calcProtonmotiveforce()

plt.subplot(411)
plt.plot(np.arange(TMAX), uin)
plt.title("Spannung innen")
plt.xlabel("Zeit")
plt.ylabel("Spannung")

plt.subplot(412)
plt.plot(np.arange(TMAX), uout)
plt.title("Spannung aussen")
plt.xlabel("Zeit")
plt.ylabel("Spannung")

plt.subplot(413)
plt.plot(np.arange(TMAX), pmftseries)
plt.title("Protonmotive force")
plt.ylabel("Zeit")
plt.ylabel("Anzahl")

plt.subplot(414)
plt.plot(np.arange(TMAX), cprotIn)
plt.title("Protonenkonzentration innen")
plt.xlabel("Zeit")
plt.ylabel("Konzentration")

plt.tight_layout()
plt.show()







