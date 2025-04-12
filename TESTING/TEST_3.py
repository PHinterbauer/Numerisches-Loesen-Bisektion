import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x**2 - 28

# Bisektion mit Tracking
def bisektion_detail(a, b, tol=0.001):
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

# Berechnung
a_vals, b_vals, c_vals = bisektion_detail(0, 28)

# Plot
x = np.linspace(0, 30, 500)
y = f(x)

plt.figure(figsize=(9, 6))
plt.plot(x, y, label='f(x) = x² - 28', color='blue')
plt.axhline(0, color='gray', linestyle='--')

# Iterationen (rote Punkte)
c_y = [f(c) for c in c_vals]
plt.plot(c_vals, c_y, 'ro', label='cₙ – Näherungen')
for i, (c, y_val) in enumerate(zip(c_vals, c_y)):
    plt.text(c, y_val + 2, f'c{i}', fontsize=9, color='darkred', ha='center')

# Vertikale Linien zu c-Werten
for c in c_vals:
    plt.vlines(c, 0, f(c), colors='red', linestyles='dotted', alpha=0.6)

# Anfangsintervall
plt.axvline(a_vals[0], color='green', linestyle='--', label='a₀ = 0')
plt.axvline(b_vals[0], color='purple', linestyle='--', label='b₀ = 28')

plt.title('Bisektionsverfahren zur Berechnung von √28')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
