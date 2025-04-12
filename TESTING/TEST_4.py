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

# Daten vorbereiten
a0, b0 = 0, 28
a_vals, b_vals, c_vals = bisektion_schritte(a0, b0)

# Plot vorbereiten
x = np.linspace(0, 30, 500)
y = f(x)

fig, ax = plt.subplots(figsize=(9, 6))
line_func, = ax.plot(x, y, label='f(x) = x² - 28', color='blue')
point_cn, = ax.plot([], [], 'ro', label='cₙ – Näherung')
vlines = []
textlabels = []

ax.axhline(0, color='gray', linestyle='--')
ax.set_xlim(0, 30)
ax.set_ylim(min(y)-5, max(y)+5)
ax.set_title('Bisektionsverfahren: Annäherung an √28')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.grid(True)

# Anfangsgrenzen
ax.axvline(a_vals[0], color='green', linestyle='--', label='a₀ = 0')
ax.axvline(b_vals[0], color='purple', linestyle='--', label='b₀ = 28')
ax.legend()

# Animationsfunktion
def update(frame):
    for v in vlines:
        v.remove()
    for t in textlabels:
        t.remove()
    vlines.clear()
    textlabels.clear()

    c_slice = c_vals[:frame+1]
    y_slice = [f(c) for c in c_slice]

    point_cn.set_data(c_slice, y_slice)

    for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
        v = ax.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
        t = ax.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
        vlines.append(v)
        textlabels.append(t)

    return line_func, point_cn, *vlines, *textlabels

ani = FuncAnimation(fig, update, frames=len(c_vals), interval=1000, repeat=False)
plt.show()
