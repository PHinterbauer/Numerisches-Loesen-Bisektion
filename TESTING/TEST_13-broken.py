import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Plotter:
    def __init__(self):
        # Attribute für die benötigten Werte
        self.a_vals = []
        self.b_vals = []
        self.c_vals = []
        self.fc_vals = []
        self.log_fc_vals = []
        self.x_dense = np.linspace(-2, 30, 500)
        self.y_dense = [0] * len(self.x_dense)  # Initialisiere mit Nullen

        # Erstelle die Figure und Subplots
        self.fig = plt.figure(figsize=(14, 8))
        self.ax1 = self.fig.add_subplot(2, 1, 1)
        self.ax2 = self.fig.add_subplot(2, 2, 3)
        self.ax3 = self.fig.add_subplot(2, 2, 4)

        # Speicher für dynamische Objekte
        self.vlines = []
        self.textlabels = []

        # Initialisiere Subplot 1
        self.line_func, = self.ax1.plot(self.x_dense, self.y_dense, label='f(x)', color='blue')
        self.point_cn, = self.ax1.plot([], [], 'ro', label='cₙ – Näherung')
        self.vline_a = self.ax1.axvline(x=0, color='green', linestyle='--', label='aₙ')
        self.vline_b = self.ax1.axvline(x=0, color='green', linestyle='--', label='bₙ')
        self.ax1.axhline(0, color='gray', linestyle='--')
        self.ax1.set_title('Bisektionsverfahren: Annäherung an √28')
        self.ax1.set_xlabel('x')
        self.ax1.set_ylabel('f(x)')
        self.ax1.grid(True)
        self.ax1.legend(loc='upper left')

        # Initialisiere Subplot 2
        self.line_log_fc, = self.ax2.plot([], [], label='|f(cₙ)|', color='green', marker='o')
        self.ax2.axhline(0, color='gray', linestyle='--')
        self.ax2.set_title('Bisektion – |f(cₙ)| pro Iteration (Log-Plot)')
        self.ax2.set_xlabel('Iteration')
        self.ax2.set_ylabel('|f(cₙ)|')
        self.ax2.set_yscale('log')
        self.ax2.grid(True)
        self.ax2.legend(loc='upper right')

        # Initialisiere Subplot 3
        self.line_fc, = self.ax3.plot([], [], label='f(cₙ)', color='green', marker='o')
        self.ax3.axhline(0, color='gray', linestyle='--')
        self.ax3.set_title('Bisektion – f(cₙ) pro Iteration (Normaler Plot)')
        self.ax3.set_xlabel('Iteration')
        self.ax3.set_ylabel('f(cₙ)')
        self.ax3.grid(True)
        self.ax3.legend(loc='upper right')

    def update_data(self, a_vals, b_vals, c_vals, fc_vals):
        """Aktualisiert die Daten der Klasse."""
        self.a_vals = a_vals
        self.b_vals = b_vals
        self.c_vals = c_vals
        self.fc_vals = np.abs(fc_vals) + 1e-10  # Konvertiere negative Werte zu positiven
        self.log_fc_vals = np.abs(fc_vals) + 1e-10
        self.y_dense = [f(x) for x in self.x_dense]
        self.line_func.set_data(self.x_dense, self.y_dense)  # Aktualisiere die Funktion f(x)

    def update(self, frame):
        """Animationsfunktion."""
        for v in self.vlines:
            v.remove()
        for t in self.textlabels:
            t.remove()
        self.vlines.clear()
        self.textlabels.clear()

        a_n, b_n = self.a_vals[frame], self.b_vals[frame]
        c_n = self.c_vals[frame]

        # Rote Punkte
        c_slice = self.c_vals[:frame+1]
        y_slice = [f(c) for c in c_slice]
        self.point_cn.set_data(c_slice, y_slice)

        # Vertikale Linien a_n, b_n
        self.vline_a.set_xdata([a_n])
        self.vline_b.set_xdata([b_n])

        # Dynamisches Zoomen in x-Richtung mit minimalem Zoom
        margin = (b_n - a_n) * 0.05
        self.ax1.set_xlim(a_n - margin, b_n + margin)

        # Vertikale Linie c_n + Text
        for i, (c, yval) in enumerate(zip(c_slice, y_slice)):
            v = self.ax1.vlines(c, 0, yval, colors='red', linestyles='dotted', alpha=0.6)
            t = self.ax1.text(c, yval + 2, f'c{i}', fontsize=9, color='darkred', ha='center')
            self.vlines.append(v)
            self.textlabels.append(t)

        # Animation für Subplot 2
        log_fc_slice = self.log_fc_vals[:frame+1]
        self.line_log_fc.set_data(range(len(log_fc_slice)), log_fc_slice)
        self.ax2.set_xlim(0, len(self.log_fc_vals))
        self.ax2.set_ylim(min(self.log_fc_vals) * 0.9, max(self.log_fc_vals) * 1.1)

        # Animation für Subplot 3
        fc_slice = self.fc_vals[:frame+1]
        self.line_fc.set_data(range(len(fc_slice)), fc_slice)
        self.ax3.set_xlim(0, len(self.fc_vals))
        self.ax3.set_ylim(min(self.fc_vals) * 1.1, max(self.fc_vals) * 1.1)

        return self.line_func, self.point_cn, self.vline_a, self.vline_b, *self.vlines, *self.textlabels

    def start_animation(self):
        """Startet die Animation."""
        ani = FuncAnimation(self.fig, self.update, frames=len(self.c_vals), interval=1000, repeat=False)
        plt.tight_layout(pad=2.0)
        plt.show()

# Bisektionsfunktion
def bisektion(f, a, b, tol, plotter):
    """Führt das Bisektionsverfahren durch und aktualisiert die Plotter-Daten."""
    a_vals, b_vals, c_vals, fc_vals = [], [], [], []
    while (b - a) / 2.0 > tol:
        c = (a + b) / 2.0
        a_vals.append(a)
        b_vals.append(b)
        c_vals.append(c)
        fc_vals.append(f(c))
        if f(c) == 0:
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
        # Aktualisiere die Daten des Plotters
        plotter.update_data(a_vals, b_vals, c_vals, fc_vals)
    return a_vals, b_vals, c_vals, fc_vals

# Beispiel für die Verwendung der Klasse
def f(x):
    return x**2 - 28

# Plotter-Instanz erstellen
plotter = Plotter()

# Bisektionsverfahren ausführen und Animation starten
bisektion(f, 0, 28, 0.001, plotter)
plotter.start_animation()