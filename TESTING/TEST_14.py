import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Plotter:
    def __init__(self):
        self.a_vals, self.b_vals, self.c_vals, self.fc_vals = [], [], [], []
        self.x_dense = np.linspace(-2, 30, 500)

        # Korrektur: plt.subplots gibt ein Array zurück
        self.fig, axes = plt.subplots(2, 2, figsize=(10, 6))
        self.ax1 = axes[0, 0]
        self.ax2 = axes[1, 0]
        self.ax3 = axes[1, 1]

        self.line_func, = self.ax1.plot([], [], label='f(x)', color='blue')
        self.point_cn, = self.ax1.plot([], [], 'ro', label='cₙ')
        self.ax1.axhline(0, color='gray', linestyle='--')
        self.ax1.legend()

        self.line_log_fc, = self.ax2.plot([], [], label='|f(cₙ)|', color='green', marker='o')
        self.ax2.set_yscale('log')
        self.ax2.legend()

        self.line_fc, = self.ax3.plot([], [], label='f(cₙ)', color='green', marker='o')
        self.ax3.legend()

    def update_data(self, a_vals, b_vals, c_vals, fc_vals):
        self.a_vals, self.b_vals, self.c_vals = a_vals, b_vals, c_vals
        self.fc_vals = np.abs(fc_vals) + 1e-10
        self.line_func.set_data(self.x_dense, [f(x) for x in self.x_dense])

    def update(self, frame):
        c_slice = self.c_vals[:frame+1]
        self.point_cn.set_data(c_slice, [f(c) for c in c_slice])
        self.line_log_fc.set_data(range(len(c_slice)), self.fc_vals[:frame+1])
        self.line_fc.set_data(range(len(c_slice)), self.fc_vals[:frame+1])
        return self.line_func, self.point_cn, self.line_log_fc, self.line_fc

    def start_animation(self):
        FuncAnimation(self.fig, self.update, frames=len(self.c_vals), interval=500, repeat=False)
        plt.show()

def bisektion(f, a, b, tol, plotter):
    a_vals, b_vals, c_vals, fc_vals = [], [], [], []
    while (b - a) / 2 > tol:
        c = (a + b) / 2
        a_vals.append(a)
        b_vals.append(b)
        c_vals.append(c)
        fc_vals.append(f(c))
        if f(a) * f(c) < 0: b = c
        else: a = c
        plotter.update_data(a_vals, b_vals, c_vals, fc_vals)
    return a_vals, b_vals, c_vals, fc_vals

def f(x): return x**2 - 28

plotter = Plotter()
bisektion(f, 0, 28, 0.001, plotter)
plotter.start_animation()