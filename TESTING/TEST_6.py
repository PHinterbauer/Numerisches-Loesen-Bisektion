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

# Startwerte
a0, b0 = 0, 28
a_vals, b_vals, c_vals = bisektion_schritte(a0, b0)

# Funktionswerte für Plot
x_dense = np.linspace(-2, 30, 500)
y_dense = f(x_dense)

# Plot-Grundgerüst
fig, ax = plt.subplots(figsize=(10, 6))
line_func, = ax.plot(x_dense, y_dense, label='f(x) = x² - 28', color='blue')
point_cn, = ax.plot([], [], 'ro', label='cₙ – Näherung')
vline_a = ax.axvline(x=0, color='green', linestyle='--', label='aₙ')
vline_b = ax.axvline(x=0, color='green', linestyle='--', label='bₙ')

ax.axhline(0, color='gray', linestyle='--')
ax.set_ylim(min(y_dense) - 5, max(y_dense) + 5)
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

    a_n, b_n = a_vals[frame], b_vals[frame]
    c_n = c_vals[frame]
    f_c_n = f(c_n)

    # Rote Punkte
    c_slice = c_vals[:frame+1]
    y_slice = [f(c) for c in c_slice]
    point_cn.set_data(c_slice, y_slice)

    # Vertikale Linien a_n, b_n
    vline_a.set_xdata([a_n])  # Hier die Liste verwenden
    vline_b.set_xdata([b_n])  # Hier die Liste verwenden

    # Dynamisches Zoomen in x-Richtung
    margin = (b_n - a_n) * 0.2
    ax.set_xlim(a_n - margin, b_n + margin)

    # Vertikale Linie c_n + Text
    for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
        v = ax.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
        t = ax.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
        vlines.append(v)
        textlabels.append(t)

    return line_func, point_cn, vline_a, vline_b, *vlines, *textlabels


ani = FuncAnimation(fig, update, frames=len(c_vals), interval=1000, repeat=False)
plt.show()
