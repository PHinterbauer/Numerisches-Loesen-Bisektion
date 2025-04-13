import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define the function
def f(x):
    return x**2 - 28

# Bisection method
def bisektion_schritte(a, b, tol=0.001):
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
    return a_vals, b_vals, c_vals, fc_vals

class Plotter:
    def __init__(self, a_vals, b_vals, c_vals, fc_vals):
        self.a_vals, self.b_vals, self.c_vals, self.fc_vals = a_vals, b_vals, c_vals, fc_vals
        self.log_fc_vals = np.abs(fc_vals) + 1e-10  # Avoid log(0)
        self.x_dense = np.linspace(-2, 30, 500)
        self.y_dense = f(self.x_dense)

        # Create figure and subplots
        self.fig = plt.figure(figsize=(14, 8))
        self._setup_subplots()
        print(self.fc_vals, self.c_vals)
        print(self.a_vals)
    def _setup_subplots(self):
        # Subplot 1: Function and bisection points
        self.ax1 = self.fig.add_subplot(2, 1, 1)
        self.line_func, = self.ax1.plot(self.x_dense, self.y_dense, label='f(x) = x² - 28', color='blue')
        self.point_cn, = self.ax1.plot([], [], 'ro', label='cₙ – Approximation')
        self.vline_a = self.ax1.axvline(0, color='green', linestyle='--', label='aₙ')
        self.vline_b = self.ax1.axvline(0, color='green', linestyle='--', label='bₙ')
        self.ax1.axhline(0, color='gray', linestyle='--')
        self.ax1.set_title('Bisection Method: Approaching √28')
        self.ax1.set_xlabel('x')
        self.ax1.set_ylabel('f(x)')
        self.ax1.grid(True)
        self.ax1.legend(loc='upper left')

        # Subplot 2: Logarithmic |f(cₙ)|
        self.ax2 = self.fig.add_subplot(2, 2, 3)
        self.line_log_fc, = self.ax2.plot([], [], 'go-', label='|f(cₙ)|')
        self.ax2.axhline(0, color='gray', linestyle='--')
        self.ax2.set_title('Bisection – |f(cₙ)| per Iteration (Log-Plot)')
        self.ax2.set_xlabel('Iteration')
        self.ax2.set_ylabel('|f(cₙ)|')
        self.ax2.set_yscale('log')
        self.ax2.grid(True)
        self.ax2.legend(loc='upper right')

        # Subplot 3: Linear f(cₙ)
        self.ax3 = self.fig.add_subplot(2, 2, 4)
        self.line_fc, = self.ax3.plot([], [], 'go-', label='f(cₙ)')
        self.ax3.axhline(0, color='gray', linestyle='--')
        self.ax3.set_title('Bisection – f(cₙ) per Iteration')
        self.ax3.set_xlabel('Iteration')
        self.ax3.set_ylabel('f(cₙ)')
        self.ax3.grid(True)
        self.ax3.legend(loc='upper right')

        # Store text objects for iteration numbers
        self.texts = []

    def update(self, frame):
        # Clear previous iteration numbers
        for text in self.texts:
            text.remove()
        self.texts.clear()

        # Update Subplot 1
        a_n, b_n, c_n = self.a_vals[frame], self.b_vals[frame], self.c_vals[frame]
        c_slice = self.c_vals[:frame+1]
        y_slice = [f(c) for c in c_slice]
        self.point_cn.set_data(c_slice, y_slice)
        self.vline_a.set_xdata([a_n])  # Wrap in a list to make it a sequence
        self.vline_b.set_xdata([b_n])  # Wrap in a list to make it a sequence
        margin = (b_n - a_n) * 0.05
        self.ax1.set_xlim(a_n - margin, b_n + margin)

        # Add iteration numbers above red points with a higher offset
        for i, (c, y) in enumerate(zip(c_slice, y_slice)):
            text = self.ax1.text(c, y + 20, f'{i}', color='red', fontsize=8, ha='center')  # Increased offset to `y + 2`
            self.texts.append(text)

        # Update Subplot 2
        log_fc_slice = self.log_fc_vals[:frame+1]
        self.line_log_fc.set_data(range(len(log_fc_slice)), log_fc_slice)
        self.ax2.set_xlim(0, len(self.log_fc_vals))
        self.ax2.set_ylim(min(self.log_fc_vals) * 0.9, max(self.log_fc_vals) * 1.1)

        # Update Subplot 3
        fc_slice = self.fc_vals[:frame+1]
        self.line_fc.set_data(range(len(fc_slice)), fc_slice)
        self.ax3.set_xlim(0, len(self.fc_vals))
        self.ax3.set_ylim(min(self.fc_vals) * 1.1, max(self.fc_vals) * 1.1)

        return self.point_cn, self.vline_a, self.vline_b, self.line_log_fc, self.line_fc

    def animate(self):
        ani = FuncAnimation(self.fig, self.update, frames=len(self.c_vals), interval=1000, repeat=False)
        plt.tight_layout()
        plt.show()

# Main execution
if __name__ == "__main__":
    a_vals, b_vals, c_vals, fc_vals = bisektion_schritte(0, 28)
    plotter = Plotter(a_vals, b_vals, c_vals, fc_vals)
    plotter.animate()