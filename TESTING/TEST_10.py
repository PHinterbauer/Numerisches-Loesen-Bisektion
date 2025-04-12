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

# Diagramm: f(cₙ) pro Iteration (Log-Plot) und f(cₙ) pro Iteration (normaler Plot)
plt.figure(figsize=(14, 4))

# Plot für f(cₙ) pro Iteration (Log-Plot)
plt.subplot(1, 2, 1)
# Werte müssen positiv sein, um sie log zu skalieren
log_fc_vals = np.abs(fc_vals)  # Absolute Werte für den Log-Plot
plt.plot(log_fc_vals, label='|f(cₙ)|', color='green', marker='o')
plt.axhline(0, color='gray', linestyle='--')
plt.title('Bisektion – |f(cₙ)| pro Iteration (Log-Plot)')
plt.xlabel('Iteration')
plt.ylabel('|f(cₙ)|')
plt.yscale('log')  # Logarithmische Skala
plt.grid(True)

# Plot für f(cₙ) pro Iteration (normaler Plot)
plt.subplot(1, 2, 2)
plt.plot(fc_vals, label='f(cₙ)', color='green', marker='o')
plt.axhline(0, color='gray', linestyle='--')
plt.title('Bisektion – f(cₙ) pro Iteration (Normaler Plot)')
plt.xlabel('Iteration')
plt.ylabel('f(cₙ)')
plt.grid(True)

plt.tight_layout()
plt.show()
