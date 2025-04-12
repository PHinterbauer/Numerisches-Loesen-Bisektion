import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x**3 - x - 2

# Bisektion mit Aufzeichnung
def bisektion_detail(a, b, tol=1e-5):
    a_vals = []
    b_vals = []
    c_vals = []
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

# Berechnungen
a_vals, b_vals, c_vals = bisektion_detail(1, 2)

# Plot vorbereiten
x = np.linspace(0, 3, 400)
y = f(x)

plt.figure(figsize=(8, 6))
plt.plot(x, y, label='f(x)', color='blue')
plt.axhline(0, color='gray', linestyle='--')

# Punkte einzeichnen: c-Werte (Nullstellenkandidaten)
plt.plot(c_vals, [f(c) for c in c_vals], 'ro', label='cₙ – Näherungen')
for c in c_vals:
    plt.vlines(c, 0, f(c), colors='red', linestyles='dotted', alpha=0.6)

# Anfangsintervall hervorheben
plt.axvline(a_vals[0], color='green', linestyle='--', label='a₀ (Start links)')
plt.axvline(b_vals[0], color='purple', linestyle='--', label='b₀ (Start rechts)')

plt.title('Bisektionsverfahren auf f(x)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
