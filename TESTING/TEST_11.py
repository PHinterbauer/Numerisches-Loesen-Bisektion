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

# Sicherstellen, dass log_fc_vals keine null- oder negativen Werte enthält
epsilon = 1e-10  # Kleiner positiver Wert
log_fc_vals = np.abs(fc_vals) + epsilon  # Konvertiere negative Werte zu positiven und füge epsilon hinzu

# Erstelle eine Figure mit mehreren Subplots
fig = plt.figure(figsize=(14, 8))  # Reduzierte Höhe der Figur

# Subplot 1: Animation des Bisektionsverfahrens (über gesamte Breite)
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

# Subplot 2: Animation für |f(cₙ)| pro Iteration
ax2 = fig.add_subplot(2, 2, 3)
line_log_fc, = ax2.plot([], [], label='|f(cₙ)|', color='green', marker='o')
ax2.axhline(0, color='gray', linestyle='--')
ax2.set_title('Bisektion – |f(cₙ)| pro Iteration (Log-Plot)')
ax2.set_xlabel('Iteration')
ax2.set_ylabel('|f(cₙ)|')
ax2.set_yscale('log')  # Logarithmische Skala
ax2.grid(True)
ax2.legend(loc='upper right')

# Subplot 3: Animation für f(cₙ) pro Iteration
ax3 = fig.add_subplot(2, 2, 4)
line_fc, = ax3.plot([], [], label='f(cₙ)', color='green', marker='o')
ax3.axhline(0, color='gray', linestyle='--')
ax3.set_title('Bisektion – f(cₙ) pro Iteration (Normaler Plot)')
ax3.set_xlabel('Iteration')
ax3.set_ylabel('f(cₙ)')
ax3.grid(True)
ax3.legend(loc='upper right')

# Animationsfunktion
def update(frame):
    # Animation für Subplot 1
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
    margin = (b_n - a_n) * 0.05  # Sehr kleine Schritte für den Zoom
    ax1.set_xlim(a_n - margin, b_n + margin)

    # Vertikale Linie c_n + Text
    for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
        v = ax1.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
        t = ax1.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
        vlines.append(v)
        textlabels.append(t)

    # Animation für Subplot 2
    log_fc_slice = log_fc_vals[:frame+1]
    line_log_fc.set_data(range(len(log_fc_slice)), log_fc_slice)
    ax2.set_xlim(0, len(log_fc_vals))
    ax2.set_ylim(min(log_fc_vals) * 0.9, max(log_fc_vals) * 1.1)

    # Animation für Subplot 3
    fc_slice = fc_vals[:frame+1]
    line_fc.set_data(range(len(fc_slice)), fc_slice)
    ax3.set_xlim(0, len(fc_vals))
    ax3.set_ylim(min(fc_vals) * 1.1, max(fc_vals) * 1.1)

    return line_func, point_cn, vline_a, vline_b, line_log_fc, line_fc, *vlines, *textlabels

ani = FuncAnimation(fig, update, frames=len(c_vals), interval=1000, repeat=False)

plt.tight_layout(pad=2.0)  # Optimiertes Layout, um White Space zu minimieren
plt.show()