import matplotlib.pyplot as plt
import copy
import numpy as np


def getDistancesAndElevationsSpace(fileName):
    distances = []
    elevations = []

    with open(fileName, 'r') as file:
        for line in file:
            values = line.split()
            wartosc1 = float(values[0])
            distances.append(wartosc1)
            wartosc2 = float(values[1])
            elevations.append(wartosc2)

    return distances, elevations


def getDistancesAndElevationsComma(fileName):
    distances = []
    elevations = []

    with open(fileName, 'r') as file:
        for line in file:
            values = line.split(',')
            wartosc1 = float(values[0])
            distances.append(wartosc1)
            wartosc2 = float(values[1])
            elevations.append(wartosc2)

    return distances, elevations


def indeksyWezlowDlaPodprzedzialow(dlugosc, ilosc_podprzedzialow):
    ilosc_wezlow = ilosc_podprzedzialow + 1
    wyniki = []
    for i in range(ilosc_wezlow):
        indeks = ((dlugosc-1)/(ilosc_wezlow-1))*i
        wyniki.append(int(indeks))
    return wyniki


def indeksyWezlow(dlugosc, ilosc_wezlow):
    wyniki = []
    for i in range(ilosc_wezlow):
        indeks = ((dlugosc - 1) / (ilosc_wezlow - 1)) * i
        wyniki.append(int(indeks))
    return wyniki


def interpolacjaPunktu(wejscie, wyjscie, x):
    wynik = 0
    for i in range(len(wejscie)):
        licznik = 1
        mianownik = 1
        for j in range(len(wejscie)):
            if j != i:
                licznik = licznik * (x - wejscie[j])
                mianownik = mianownik * (wejscie[i] - wejscie[j])
        skladnik = (licznik/mianownik) * wyjscie[i]
        wynik = wynik + skladnik
    return wynik


def wygenerujWartosciInterpolacji(wejscie, wyjscie, indeksy):
    tabela = []
    for i in range(len(wejscie)):
        interpolowany = interpolacjaPunktu([wejscie[k] for k in indeksy], [wyjscie[k] for k in indeksy], wejscie[i])
        tabela.append(interpolowany)
    return tabela


def wygenerujWartosciInterpolacjiMod(wejscie, wyjscie, indeksyMod):
    tabela = []
    for i in range(len(wejscie)):
        interpolowany = interpolacjaPunktu([wejscie[k] for k in indeksyMod], [wyjscie[k] for k in indeksyMod], wejscie[i])
        tabela.append(interpolowany)
    return tabela


def zwrocInterpolacje(wejscie, wyjscie, ilosc_wezlow):
    indeksy = indeksyWezlow(len(wejscie), ilosc_wezlow)
    # wynik = interpolacjaPunktu([wejscie[k] for k in indeksy], [wyjscie[k] for k in indeksy], wejscie[20])
    interpolowane = wygenerujWartosciInterpolacji(wejscie, wyjscie, indeksy)
    return interpolowane

def zwrocInterpolacjeZmodyfkowane(wejscie, wyjscie, ilosc_wezlow):
    indeksyMod = [0, 10, 20, 60, 100, 140, 180, 220, 240, 260, 300, 340, 380, 420, 460, 470, 480]
    # wynik = interpolacjaPunktu([wejscie[k] for k in indeksy], [wyjscie[k] for k in indeksy], wejscie[20])
    interpolowane = wygenerujWartosciInterpolacjiMod(wejscie, wyjscie, indeksyMod)
    return interpolowane


def wyzanczMacierzInterpolacji(punkty, odl):
    macierz_glowna = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [1, odl, odl ** 2, odl ** 3, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1, odl, odl ** 2, odl ** 3],
    [0, 1, 2*odl, 3*(odl**2), 0, -1, 0, 0],
    [0, 0, 2, 6*odl, 0, 0, -2, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 2, 6*odl]
    ]

    macierz_glowna_array = np.array(macierz_glowna)

    macierz_wynikowa = [punkty[0], punkty[1], punkty[1], punkty[2], 0, 0, 0, 0]
    macierz_wynikowa_array = np.array(macierz_wynikowa)

    macierz_wspolczynikow_array = np.linalg.solve(macierz_glowna_array, macierz_wynikowa_array)
    macierz_wspolczynikow = macierz_wspolczynikow_array.tolist()

    return macierz_wspolczynikow


def wyznaczWspolczynniki(punkty, h, ilosc):
    tab = [[0 for a in range(ilosc*4)] for b in range(ilosc*4)]

    for i in range(ilosc*2):
        if i%2 == 0:
            tab[i][2*i] = 1
        else:
            tab[i][2 * (i - 1) + 0] = 1
            tab[i][2 * (i - 1) + 1] = h
            tab[i][2 * (i - 1) + 2] = h ** 2
            tab[i][2 * (i - 1) + 3] = h ** 3

    for i in range(ilosc*2, ilosc*2+ilosc-1):
        wzgledny = i - ilosc*2
        tab[i][wzgledny * 4 + 1] = 1
        tab[i][wzgledny * 4 + 2] = 2*h
        tab[i][wzgledny * 4 + 3] = 3*(h**2)
        tab[i][wzgledny * 4 + 5] = -1

    for i in range(ilosc*2+ilosc-1, ilosc*4-2):
        wzgledny = i - (ilosc*2+ilosc-1)
        tab[i][wzgledny * 4 + 2] = 2
        tab[i][wzgledny * 4 + 3] = 6 * h
        tab[i][wzgledny * 4 + 6] = -2

    tab[ilosc * 4 - 2][2] = 1
    tab[ilosc * 4 - 1][ilosc * 4 - 2] = 2
    tab[ilosc * 4 - 1][ilosc * 4 - 1] = 6*h

    tab_array = np.array(tab)

    wynikowa = [0 for c in range(ilosc*4)]
    licznik = 0
    for i in range(ilosc*2):
        if i%2 == 0:
            wynikowa[i] = punkty[licznik]
            licznik = licznik + 1
        else:
            wynikowa[i] = punkty[licznik]

    wynikowa_array = np.array(wynikowa)

    wspolczynnikow_array = np.linalg.solve(tab_array, wynikowa_array)
    wspolczynnikow = wspolczynnikow_array.tolist()

    return wspolczynnikow


def okreslWspolczynniki(wejscie, wyjscie, ilosc_przedzialow):
    indeksy = indeksyWezlowDlaPodprzedzialow(len(wejscie), ilosc_przedzialow)
    punkty = [wyjscie[k] for k in indeksy]
    odl = wejscie[indeksy[1]] - wejscie[indeksy[0]]
    ilosc = ilosc_przedzialow
    wsp = wyznaczWspolczynniki(punkty, odl, ilosc)
    return wsp


def okreslFunkcjeInterpolowana(wejscie, wyjscie, ilosc_przedzialow):
    wsp = okreslWspolczynniki(wejscie, wyjscie, ilosc_przedzialow)
    indeksy = indeksyWezlowDlaPodprzedzialow(len(wejscie), ilosc_przedzialow)
    wyniki = []
    for i in range(len(wejscie)):
        iloraz = ((i*(ilosc_przedzialow))/(len(wejscie)))
        glowny = int(iloraz)
        a, b, c, d = wsp[4 * glowny], wsp[4 * glowny + 1], wsp[4 * glowny + 2], wsp[4 * glowny + 3]
        roznica = (wejscie[i] - wejscie[indeksy[glowny]])
        posredni = a + b * roznica + c * roznica * roznica + d * roznica * roznica * roznica
        wyniki.append(posredni)
    return wyniki


dystanse, wysokosci = getDistancesAndElevationsSpace('chelm.txt')
print(dystanse)
print(wysokosci)
print(len(dystanse))
print(len(wysokosci))
maksimum = max(wysokosci)
minimum = min(wysokosci)
print("Różnica wysokości to", maksimum-minimum)


plt.plot(dystanse, wysokosci)
plt.title('Wykres wysokości trasy Chełm')
plt.xlabel('Dystans w metrach')
plt.ylabel('Wysokość w metrach')
plt.savefig('ChelmOryginalny.png')
plt.show()


'''
indeksy = indeksyWezlow(len(dystanse), 8)
print(indeksy)
wynik = interpolacjaPunktu([dystanse[k] for k in indeksy], [wysokosci[k] for k in indeksy], dystanse[20])
print(wynik)
print(wysokosci[20])
interpolowane = wygenerujWartosciInterpolacji(dystanse, wysokosci, indeksy)
print(interpolowane)
'''

'''
# WAZNE
iloscWezlow = 17

interpolowane = zwrocInterpolacje(dystanse, wysokosci, iloscWezlow)
indeksy = indeksyWezlow(len(dystanse), iloscWezlow)

plt.plot(dystanse, wysokosci)
plt.plot(dystanse, interpolowane)
plt.scatter([dystanse[k] for k in indeksy], [wysokosci[k] for k in indeksy], label='Interpolated Points', color='red')
plt.title("Interpolacja trasy Chełm metodą Lagrange'a dla {} węzłów".format(iloscWezlow))
plt.xlabel('Dystans w metrach')
plt.ylabel('Wysokość w metrach')
plt.legend(["Funkcja interpolowana", "Funkcja interpolacyjna"])
plt.savefig("ChelmLagrange{}Wezly.png".format(iloscWezlow))
plt.show()
'''

# WAZNE
iloscWezlow = 17

interpolowane = zwrocInterpolacjeZmodyfkowane(dystanse, wysokosci, iloscWezlow)
indeksyMod = [0, 10, 20, 60, 100, 140, 180, 220, 240, 260, 300, 340, 380, 420, 460, 470, 480]

plt.plot(dystanse, wysokosci)
plt.plot(dystanse, interpolowane)
plt.scatter([dystanse[k] for k in indeksyMod], [wysokosci[k] for k in indeksyMod], label='Interpolated Points', color='red')
plt.title("Interpolacja trasy Chełm zmodyfikowaną metodą Lagrange'a dla {} węzłów".format(iloscWezlow))
plt.xlabel('Dystans w metrach')
plt.ylabel('Wysokość w metrach')
plt.legend(["Funkcja interpolowana", "Funkcja interpolacyjna"])
plt.savefig("ChelmZmodyfikowanyLagrange{}Wezly.png".format(iloscWezlow))
plt.show()


'''
punkty = [6, -2, 4]
odleglosc = 2
wspolczynniki = wyzanczMacierzInterpolacji(punkty, odleglosc)
print(wspolczynniki)

punkty = [-2, 4, 9]
odleglosc = 2
wspolczynniki = wyzanczMacierzInterpolacji(punkty, odleglosc)
print(wspolczynniki)
'''

'''
punkty = [21, 24, 24, 18, 16]
odleglosc = 1
wspolczynniki = wyznaczWspolczynniki(punkty, odleglosc, len(punkty)-1)
print(wspolczynniki)

wyniki = []
for i in np.arange(0, 4, 0.001):
    if i <= 0+1:
        posredni = wspolczynniki[0] + wspolczynniki[1] * (i-0) + wspolczynniki[2] * (i-0) * (i-0) + wspolczynniki[3] * (i-0) * (i-0) * (i-0)
        wyniki.append(posredni)
    elif i < 2:
        posredni = wspolczynniki[4] + wspolczynniki[5] * (i - 1) + wspolczynniki[6] * (i - 1) * (i - 1) + wspolczynniki[7] * (i - 1) * (i - 1) * (i - 1)
        wyniki.append(posredni)
    elif i < 3:
        posredni = wspolczynniki[8] + wspolczynniki[9] * (i - 2) + wspolczynniki[10] * (i - 2) * (i - 2) + wspolczynniki[11] * (i - 2) * (i - 2) * (i - 2)
        wyniki.append(posredni)
    elif i < 4:
        posredni = wspolczynniki[12] + wspolczynniki[13] * (i - 3) + wspolczynniki[14] * (i - 3) * (i - 3) + wspolczynniki[15] * (i - 3) * (i - 3) * (i - 3)
        wyniki.append(posredni)

plt.plot(np.arange(0, 4, 0.001), wyniki)
plt.show()
'''

'''
# WAZNE
iloscPodrzedzialow = 96

wspolczynniki = okreslWspolczynniki(dystanse, wysokosci, iloscPodrzedzialow)
print(wspolczynniki)

daneInterpolowane = okreslFunkcjeInterpolowana(dystanse, wysokosci, iloscPodrzedzialow)
print(daneInterpolowane)

indeksy = indeksyWezlowDlaPodprzedzialow(len(dystanse), iloscPodrzedzialow)

plt.plot(dystanse, wysokosci)
plt.plot(dystanse, daneInterpolowane)
plt.scatter([dystanse[k] for k in indeksy], [wysokosci[k] for k in indeksy], label='Interpolated Points', color='red')
plt.title("Interpolacja trasy Ulm-Lugano metodą splajnów sześciennych dla {} węzłów".format(iloscPodrzedzialow+1))
plt.xlabel('Dystans w metrach')
plt.ylabel('Wysokość w metrach')
plt.legend(["Funkcja interpolowana", "Funkcja interpolacyjna"])
plt.savefig("UlmLuganoSplajnySzescienne{}Wezly.png".format(iloscPodrzedzialow+1))
plt.show()
'''
