import math

# Polynomdefinition
def f(x):
    return 2 * x + x**2 + 3 * x**3 - x**4

# Intervall finden mit Vorzeichenwechsel um x0
def finde_intervall(f, x0, max_versuche=10000000000000000000, startweite=1e-4, schrittweite=1e-4):
    a = x0 - startweite
    b = x0 + startweite
    for _ in range(max_versuche):
        if f(a) * f(b) < 0:
            return a, b
        a -= schrittweite
        b += schrittweite
    raise ValueError("Kein gültiges Intervall mit Vorzeichenwechsel gefunden.")

# Bisektion mit Ausgabe + Rückgabe des letzten sinnvollen Intervalls
def bisektion_mit_druck(f, a, b, tol=1e-10, max_iter=500000000000):
    print("\n--- Iterationen der Bisektionsmethode ---")
    print(f"{'Nr.':>3} | {'a':>12} | {'b':>12} | {'c':>12} | {'f(c)':>12}")
    print("-" * 60)
    for i in range(max_iter):
        c = (a + b) / 2
        fc = f(c)
        print(f"{i:>3} | {a:12.8f} | {b:12.8f} | {c:12.8f} | {fc:12.8f}")
        if abs(fc) < tol or abs(b - a) < tol:
            # Endintervall bestimmen
            if f(a) * fc < 0:
                return c, i + 1, a, c
            else:
                return c, i + 1, c, b
        if f(a) * fc < 0:
            b = c
        else:
            a = c
    # Falls maximale Iterationen erreicht
    if f(a) * f((a + b) / 2) < 0:
        return (a + b) / 2, max_iter, a, (a + b) / 2
    else:
        return (a + b) / 2, max_iter, (a + b) / 2, b

# Hauptprogramm
x0 = 3.4567
a, b = finde_intervall(f, x0)
epsilon = abs(b - a)
potenz = round(math.log10(epsilon / abs(x0)))

print(f"Gefundenes Intervall: a = {a}, b = {b}")
print(f"f(a) = {f(a):.6f}, f(b) = {f(b):.6f}")
print(f"Intervallweite (epsilon): {epsilon:.2e}")
print(f"Genauigkeit: 10^{potenz}")

wurzel, anzahl_iterationen, intervall1, intervall2 = bisektion_mit_druck(f, a, b)

print(f"\nGefundene Nullstelle: x = {wurzel:.10f} nach {anzahl_iterationen} Iterationen")
print(f"Letztes Intervall: [{intervall1:.10f}, {intervall2:.10f}]")
