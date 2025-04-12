import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Funktion definieren
def f(x):
    return x**2 - 28

# Bisektionsverfahren mit Zwischenspeicherung
def bisektion_schritte(a, b, tol=0.001):
    a_vals, b_vals, c_vals, fc_vals, intervals = [], [], [], [], []
    while (b - a) / 2.0 > tol:
        c = (a + b) / 2.0
        a_vals.append(a)
        b_vals.append(b)
        c_vals.append(c)
        fc_vals.append(f(c))
        intervals.append(b - a)
        if f(c) == 0:
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return a_vals, b_vals, c_vals, fc_vals, intervals

# Startwerte
a0, b0 = 0, 28
a_vals, b_vals, c_vals, fc_vals, intervals = bisektion_schritte(a0, b0)

# Funktionswerte für Plot
x_dense = np.linspace(-2, 30, 500)
y_dense = f(x_dense)

# Erstelle eine Figure mit mehreren Subplots
fig = plt.figure(figsize=(14, 8))

# Subplot 1: Animation des Bisektionsverfahrens
ax1 = fig.add_subplot(2, 1, 1)
line_func, = ax1.plot(x_dense, y_dense, label='f(x) = x² - 28', color='blue')
point_cn, = ax1.plot([], [], 'ro', label='cₙ – Näherung')
vline_a = ax1.axvline(x=0, color='green', linestyle='--', label='aₙ')
vline_b = ax1.axvline(x=0, color='green', linestyle='--', label='bₙ')

ax1.axhline(0, color='gray', linestyle='--')
ax1.set_ylim(min(y_dense) - 5, max(y_dense) + 5)
ax1.set_title('Bisektionsverfahren: Annäherung an √28')
ax1.set_xlabel('x')
ax1.set_ylabel('f(x)')
ax1.grid(True)
ax1.legend(loc='upper left')

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
    vline_a.set_xdata([a_n])
    vline_b.set_xdata([b_n])

    # Dynamisches Zoomen in x-Richtung mit minimalem Zoom
    margin = (b_n - a_n) * 0.01
    if margin < 0.01:
        margin = 0.005

    ax1.set_xlim(a_n - margin, b_n + margin)

    # Vertikale Linie c_n + Text
    for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
        v = ax1.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
        t = ax1.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
        vlines.append(v)
        textlabels.append(t)

    return line_func, point_cn, vline_a, vline_b, *vlines, *textlabels

ani = FuncAnimation(fig, update, frames=len(c_vals), interval=1000, repeat=False)

# Subplot 2: Diagramme für f(cₙ) pro Iteration
ax2 = fig.add_subplot(2, 2, 3)
log_fc_vals = np.abs(fc_vals)  # Absolute Werte für den Log-Plot
ax2.plot(log_fc_vals, label='|f(cₙ)|', color='green', marker='o')
ax2.axhline(0, color='gray', linestyle='--')
ax2.set_title('Bisektion – |f(cₙ)| pro Iteration (Log-Plot)')
ax2.set_xlabel('Iteration')
ax2.set_ylabel('|f(cₙ)|')
ax2.set_yscale('log')  # Logarithmische Skala
ax2.grid(True)

ax3 = fig.add_subplot(2, 2, 4)
ax3.plot(fc_vals, label='f(cₙ)', color='green', marker='o')
ax3.axhline(0, color='gray', linestyle='--')
ax3.set_title('Bisektion – f(cₙ) pro Iteration (Normaler Plot)')
ax3.set_xlabel('Iteration')
ax3.set_ylabel('f(cₙ)')
ax3.grid(True)

plt.tight_layout()
plt.show()