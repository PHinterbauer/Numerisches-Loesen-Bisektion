import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x**3 - x - 2

def bisektion_detail(a, b, tol=1e-5):
    a_vals = []
    b_vals = []
    c_vals = []
    fc_vals = []
    intervals = []

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

# Berechnung
a_vals, b_vals, c_vals, fc_vals, intervals = bisektion_detail(1, 2)

# Diagramm 1: Intervallgrenzen a_n und b_n
plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 1)
plt.plot(a_vals, label='aₙ', marker='o')
plt.plot(b_vals, label='bₙ', marker='x')
plt.title('Bisektion – Intervallgrenzen pro Iteration')
plt.xlabel('Iteration')
plt.ylabel('Wert')
plt.legend()
plt.grid(True)

# Diagramm 2: Funktionswert f(cₙ) pro Iteration
plt.subplot(1, 3, 2)
plt.plot(fc_vals, label='f(cₙ)', color='green', marker='o')
plt.axhline(0, color='gray', linestyle='--')
plt.title('Bisektion – f(cₙ) pro Iteration')
plt.xlabel('Iteration')
plt.ylabel('f(cₙ)')
plt.grid(True)

# Diagramm 3: Intervallbreite bₙ - aₙ
plt.subplot(1, 3, 3)
plt.plot(intervals, label='Intervallbreite', color='purple', marker='s')
plt.title('Bisektion – Intervallbreite pro Iteration')
plt.xlabel('Iteration')
plt.ylabel('bₙ - aₙ')
plt.grid(True)

plt.tight_layout()
plt.show()
