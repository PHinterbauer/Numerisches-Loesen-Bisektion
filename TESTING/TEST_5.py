import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Funktion definieren
def f(x):
    return x**2 - 28

# Bisektionsverfahren mit Zwischenspeicherung
def bisektion_schritte(a, b, tol=0.001):
    a_vals, b_vals, c_vals = [], [], []
    while (b - a) / 2.0 > tol:
        c = (a + b) / 2.0
        a_vals.append(a)
        b_vals.append(b)
        c_vals.append(c)
        if f(c) == 0:
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return a_vals, b_vals, c_vals

# Intervall und Berechnung
a0, b0 = 0, 28
a_vals, b_vals, c_vals = bisektion_schritte(a0, b0)

# x-Werte für Funktionsplot
x = np.linspace(-2, 30, 500)
y = f(x)

# Plot-Grundgerüst
fig, ax = plt.subplots(figsize=(10, 6))
line_func, = ax.plot(x, y, label='f(x) = x² - 28', color='blue')
point_cn, = ax.plot([], [], 'ro', label='cₙ – Näherung')
interval_line, = ax.plot([], [], 'g-', lw=2, label='Intervall [aₙ, bₙ]')

ax.axhline(0, color='gray', linestyle='--')
ax.set_xlim(-2, 30)
ax.set_ylim(min(y)-5, max(y)+5)
ax.set_title('Bisektionsverfahren: Annäherung an √28')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.grid(True)
ax.legend(loc='upper left')

# Speicher für dynamische Objekte
vlines = []
textlabels = []

# Animationsfunktion
def update(frame):
    for v in vlines:
        v.remove()
    for t in textlabels:
        t.remove()
    vlines.clear()
    textlabels.clear()

    # Aktuelles Intervall
    a_n, b_n = a_vals[frame], b_vals[frame]
    c_n = c_vals[frame]
    f_c_n = f(c_n)

    # Intervall-Linie zeichnen
    interval_line.set_data([a_n, b_n], [0, 0])

    # Rote Punkte + Linien
    c_slice = c_vals[:frame+1]
    y_slice = [f(c) for c in c_slice]
    point_cn.set_data(c_slice, y_slice)

    for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
        v = ax.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
        t = ax.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
        vlines.append(v)
        textlabels.append(t)

    return line_func, point_cn, interval_line, *vlines, *textlabels

ani = FuncAnimation(fig, update, frames=len(c_vals), interval=1000, repeat=False)
plt.show()
