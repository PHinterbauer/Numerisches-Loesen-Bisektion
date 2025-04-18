import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Base():
    def __init__(self, 
                formula: str,
                accuracy: float,
                a: float = 0.0, 
                b: float = 0.0, 
                n: float = 0.0, 
                sep_length: int = 100, 
                methods: list = ["bisektion", "newton_raphson", "regula_falsi"], 
                enable_plot: bool = False, 
                enable_zoom: bool = False):
        self.formula = formula
        self.accuracy = accuracy
        self.a = a
        self.b = b
        self.n = n
        self.sep_length = sep_length
        self.methods = methods
        self.enable_plot = enable_plot
        self.enable_zoom = enable_zoom

    def separator(self):
        print("-" * self.sep_length)

    def get_plotter_zoom(self):
        flag = True
        while flag:
            choice = input("Do you want to activate automatic zooming on the main plot? [y/n]: ").lower()
            self.separator()
            if choice in ["y", "yes", "n", "no"]:
                if choice in ["y", "yes"]:
                    self.enable_zoom = True
                    flag = False
                else:
                    self.enable_zoom = False
                    flag = False
            else:
                print("Invalid selection! Please choose [y/n].")
                continue

    def get_plotter(self):
        flag = True
        while flag:
            choice = input("Do you want to plot the outcome? [y/n]: ").lower()
            self.separator()
            if choice in ["y", "yes", "n", "no"]:
                if choice in ["y", "yes"]:
                    self.get_plotter_zoom()
                    self.enable_plot = True
                    flag = False
                else:
                    self.enable_plot = False
                    flag = False
            else:
                print("Invalid selection! Please choose [y/n].")
                continue

    def start_method(self):
        self.separator()
        print("Select the method:\n1. Bisection\n2. Newton-Raphson\n3. Regula Falsi")
        flag = True
        while flag:
            try:
                choice = int(input("Enter the number of the method: "))
                self.separator()
                if choice in [1, 2, 3]:
                    self.get_plotter()
                    eval(self.methods[choice - 1] + ".main_loop()")
                    if self.enable_plot:
                        # plotter.update_data(eval(self.methods[choice -1]).a, eval(self.methods[choice -1]).b, eval(self.methods[choice -1]).c, eval(self.methods[choice -1]).fc, self.n)
                        plotter.plot()
                    flag = False
                else:
                    print("Invalid selection! Please choose 1, 2, or 3.")
                    continue
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_n(self):
        flag = True
        while flag:
            try:
                self.n = float(input(f'Enter n for: '))
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_interval_a(self):
        flag = True
        while flag:
            try:
                self.a = float(input(f'Enter a for: '))
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_interval_b(self):
        flag = True
        while flag:
            try:
                self.b = float(input(f'Enter b for: '))
                self.separator()
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

class BaseCalculations(Base):
    def __init__(self,
            formula: str,
            accuracy: float,
            a: float = 0.0, 
            b: float = 0.0, 
            n: float = 0.0, 
            sep_length: int = 100, 
            methods: list = ["bisektion", "newton_raphson", "regula_falsi"], 
            enable_plot: bool = False, 
            enable_zoom: bool = False, 
            c: float = 0.0, 
            fa: float = 0.0, 
            fb: float = 0.0, 
            fc: float = 0.0, 
            control_value: float = 0.0,
            result: dict = {},
            control: bool = False):
        super().__init__(formula, accuracy, a, b, n, sep_length, methods, enable_plot, enable_zoom)
        self.c = c
        self.fa = fa
        self.fb = fb
        self.fc = fc
        self.control_value = control_value
        self.result = result
        self.control = control

    def calculate_fa(self):
        self.fa = eval(self.formula.format(x=self.a, n=self.n))

    def calculate_fb(self):
        self.fb = eval(self.formula.format(x=self.b, n=self.n))

    def calculate_c(self):
        raise NotImplementedError("Die Methode 'calculate_c' muss in der Unterklasse überschrieben werden.")

    def calculate_fc(self):
        self.fc = eval(self.formula.format(x=self.c, n=self.n))

    def check_control(self):
        self.control_value = self.fa * self.fb
        self.control = self.control_value < 0

    def switch_interval(self):
        if self.fa * self.fc < 0:
            self.b = self.c
        else:
            self.a = self.c

    def calc_separator(self, text_1: str, header: str, text_2: str = "", text_3: str = ""):
        total_length = max(len(text_1), len(text_2), len(text_3))
        side_length = (total_length - len(header)) // 2

        left_side_length = side_length
        right_side_length = total_length - len(header) - left_side_length

        return left_side_length, right_side_length, total_length

    def main_loop(self):
        iteration = 1
        self.result = {}
        self.get_n()
        self.get_interval_a()
        self.get_interval_b()
        self.calculate_fa()
        self.calculate_fb()
        self.calculate_c()
        self.calculate_fc()
        self.check_control()

        if self.enable_plot:
            plotter.update_data(self.a, self.b, self.c, self.fc, self.n)

        left_side_length, right_side_length, total_length = self.calc_separator(
            text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
            text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
            text_3=f'control: {self.control}, control value: {self.control_value}',
            header=f' Iteration: {iteration} '
        )
        print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
        print(f'a: {self.a}, b: {self.b}, c: {self.c}')
        print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
        print(f'control: {self.control}, control value: {self.control_value}')
        print("=" * total_length)

        self.switch_interval()

        while np.floor(abs(self.fc) * 10 ** (abs(int(np.log10(self.accuracy))) - 1)) != 0 and self.control:
            iteration += 1
            self.calculate_fa()
            self.calculate_fb()
            self.calculate_c()
            self.calculate_fc()
            self.check_control()

            if self.enable_plot:
                plotter.update_data(self.a, self.b, self.c, self.fc, self.n)

            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
                text_3=f'control: {self.control}, control value: {self.control_value}',
                header=f' Iteration: {iteration} '
            )
            print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
            print(f'a: {self.a}, b: {self.b}, c: {self.c}')
            print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
            print(f'control: {self.control}, control value: {self.control_value}')
            print("=" * total_length)

            self.switch_interval()

            if self.fc == 0:
                self.c = self.c
                break
            elif self.fa == 0:
                self.c = self.a
                break
            elif self.fb == 0:
                self.c = self.b
                break

        if not self.control:
            print("No solution found!")
        else:
            self.result["c"] = self.c

        if self.result:
            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                header=f' Result '
            )
            print("=" * left_side_length, "Result", "=" * right_side_length)
            for key, value in self.result.items():
                result_text = f'{key}: {value}'
                left_padding, right_padding, _ = self.calc_separator(text_1="=" * total_length, header=result_text)
                print(" " * left_padding + result_text + " " * right_padding)
            print("=" * total_length)
        else:
            print("No solution found!")

class Bisektion(BaseCalculations):
    def __init__(self,
            formula,
            accuracy,
            a = 0.0, 
            b = 0.0, 
            n = 0.0, 
            sep_length = 100, 
            methods = ["bisektion", "newton_raphson", "regula_falsi"],
            enable_plot = False, 
            enable_zoom = False, 
            c = 0.0, 
            fa = 0.0, 
            fb = 0.0, 
            fc = 0.0, 
            control_value = 0.0,
            result = {},
            control = False):
        super().__init__(formula, accuracy, a, b, n, sep_length, methods, enable_plot, enable_zoom, c, fa, fb, fc, control_value, result, control)

    def calculate_c(self):
        self.c = (self.a + self.b) / 2

class NewtonRaphson(BaseCalculations):
    def __init__(self,
            formula: str,
            accuracy: float,
            formula_derivative: str,
            a: float = 0.0, 
            b: float = 0.0, 
            n: float = 0.0, 
            sep_length: int = 100, 
            methods: list = ["bisektion", "newton_raphson", "regula_falsi"],
            enable_plot: bool = False, 
            enable_zoom: bool = False, 
            c: float = 0.0, 
            fa: float = 0.0, 
            fb: float = 0.0, 
            fc: float = 0.0, 
            control_value: float = 0.0,
            result: dict = {},
            control: bool = False,
            c_derivative: float = 0.0):
        super().__init__(formula, accuracy, a, b, n, sep_length, methods, enable_plot, enable_zoom, c, fa, fb, fc, control_value, result, control)
        self.formula_derivative = formula_derivative
        self.c_derivative = c_derivative

    def calculate_c(self):
        self.derivative = eval(self.formula_derivative.format(x=self.c))
        self.c = self.c - self.fc / self.derivative

class RegulaFalsi(BaseCalculations):
    def __init__(self,
            formula,
            accuracy,
            a = 0.0, 
            b = 0.0, 
            n = 0.0, 
            sep_length = 100, 
            methods = ["bisektion", "newton_raphson", "regula_falsi"],
            enable_plot = False, 
            enable_zoom = False, 
            c = 0.0, 
            fa = 0.0, 
            fb = 0.0, 
            fc = 0.0, 
            control_value = 0.0,
            result = {},
            control = False):
        super().__init__(formula, accuracy, a, b, n, sep_length, methods, enable_plot, enable_zoom, c, fa, fb, fc, control_value, result, control)

    def calculate_c(self):
        self.c = (self.a * self.fb - self.b * self.fa) / (self.fb - self.fa)

class Plotter(Base):
    def __init__(self,
            formula: str,
            accuracy: float,
            a: float = 0.0, 
            b: float = 0.0, 
            n: float = 0.0, 
            sep_length: int = 100, 
            methods: list = ["bisektion", "newton_raphson", "regula_falsi"],
            enable_plot: bool = False, 
            enable_zoom: bool = False,
            a_list: list = [],
            b_list: list = [],
            c_list: list = [],
            fc_list: list = [],
            c_points_text_str_list: list = []):
        super().__init__(formula, accuracy, a, b, n, sep_length, methods, enable_plot, enable_zoom)
        self.a_list = a_list
        self.b_list = b_list
        self.c_list = c_list
        self.fc_list = fc_list
        self.c_points_text_str_list = c_points_text_str_list

    def initialize_plot(self, a, b, n):
        self.fig = plt.figure(figsize=(14, 8))
        self.fig.canvas.manager.set_window_title("Nullstellenberechnung durch Iterative-Verfahren")

        self.x_func_line = np.linspace(a, b, 500)
        self.y_func_line = [eval(self.formula.format(x=x, n=n)) for x in self.x_func_line]

        self.approach_to_root = self.fig.add_subplot(2, 1, 1)
        self.func_line_list, = self.approach_to_root.plot(self.x_func_line, self.y_func_line, label=self.formula, color="blue")
        self.c_point_list, = self.approach_to_root.plot([], [], "ro", label="Geschätzte Nullstelle")
        self.a_vertical_line = self.approach_to_root.axvline(0, color="green", linestyle="--", label="a")
        self.b_vertical_line = self.approach_to_root.axvline(0, color="purple", linestyle="--", label="b")
        self.approach_to_root.axhline(0, color="gray", linestyle="--")
        self.approach_to_root.set_title("Annäherung an Nullstelle")
        self.approach_to_root.set_xlabel("c")
        self.approach_to_root.set_ylabel("f(c)")
        self.approach_to_root.grid(True)
        self.approach_to_root.legend(loc="upper left")

        self.fc_per_iter_logplot = self.fig.add_subplot(2, 2, 3)
        self.fc_per_iter_logplot_line_list, = self.fc_per_iter_logplot.plot([0], [1e-100], "go-", label="|f(c)|")
        self.fc_per_iter_logplot.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter_logplot.set_title("|f(c)| pro Iteration (Logarithmisch)")
        self.fc_per_iter_logplot.set_xlabel("Iteration")
        self.fc_per_iter_logplot.set_ylabel("|f(c)|")
        self.fc_per_iter_logplot.set_yscale("log")
        self.fc_per_iter_logplot.grid(True)
        self.fc_per_iter_logplot.legend(loc="upper right")

        self.fc_per_iter = self.fig.add_subplot(2, 2, 4)
        self.fc_per_iter_line_list, = self.fc_per_iter.plot([], [], "go-", label="f(c)")
        self.fc_per_iter.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter.set_title("f(c) pro Iteration")
        self.fc_per_iter.set_xlabel("Iteration")
        self.fc_per_iter.set_ylabel("f(c)")
        self.fc_per_iter.grid(True)
        self.fc_per_iter.legend(loc="upper right")

    def update_data(self, a, b, c, fc, n):
        self.a = a
        self.b = b
        self.n = n
        self.a_list.append(self.a)
        self.b_list.append(self.b)
        self.c_list.append(c)
        self.fc_list.append(fc)

    def update_plot(self, frame):
        for text_str in self.c_points_text_str_list:
            text_str.remove()
        self.c_points_text_str_list.clear()

        a_iter = self.a_list[frame]
        b_iter = self.b_list[frame]
        c_point_x_list = self.c_list[:frame + 1]
        c_point_y_list = self.fc_list[:frame + 1]

        self.c_point_list.set_data(c_point_x_list, c_point_y_list)
        self.a_vertical_line.set_xdata([a_iter])
        self.b_vertical_line.set_xdata([b_iter])

        if self.enable_zoom:
            approach_to_root_min_width = abs(a_iter - b_iter)
            approach_to_root_zoom = approach_to_root_min_width * 0.05
            self.approach_to_root.set_xlim(a_iter - approach_to_root_zoom, b_iter + approach_to_root_zoom)

        c_point_text_y_min, c_point_text_y_max = self.approach_to_root.get_ylim()
        c_point_text_y_offset = (c_point_text_y_max - c_point_text_y_min) * 0.05
        for iter in range(len(c_point_x_list)):
            c_point_text_x_iter = c_point_x_list[iter]
            c_point_text_y_iter = c_point_y_list[iter]
            c_point_text_str_iter = self.approach_to_root.text(
                c_point_text_x_iter, c_point_text_y_iter + c_point_text_y_offset, str(iter + 1), color="red", fontsize=8, ha="center"
            )
            self.c_points_text_str_list.append(c_point_text_str_iter)

        log_fc_list = np.abs(self.fc_list[:frame + 1]) + 1e-10
        self.fc_per_iter_logplot_line_list.set_data(range(1, frame + 2), log_fc_list)
        self.fc_per_iter_logplot.set_xlim(1, len(self.c_list) + 1)
        self.fc_per_iter_logplot.set_ylim(min(log_fc_list) * 0.9, max(log_fc_list) * 1.1)

        self.fc_per_iter_line_list.set_data(range(1, frame + 2), self.fc_list[:frame + 1])
        self.fc_per_iter.set_xlim(1, len(self.c_list) + 1)
        self.fc_per_iter.set_ylim(min(self.fc_list) * 1.5, max(self.fc_list) * 1.1)

        return self.c_point_list, self.a_vertical_line, self.b_vertical_line, self.fc_per_iter_logplot_line_list, self.fc_per_iter_line_list

    def plot(self):
        self.initialize_plot(self.a, self.b, self.n)
        anim1 = animation.FuncAnimation(self.fig, self.update_plot, frames=len(self.c_list), interval=1000, repeat=False)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    formula = "{x}**2 - {n}"
    formula_derivative = "2 * {x}"
    accuracy = 0.001

    base = Base(formula, accuracy)
    bisektion = Bisektion(formula, accuracy)
    newton_raphson = NewtonRaphson(formula, accuracy, formula_derivative)
    regula_falsi = RegulaFalsi(formula, accuracy)
    plotter = Plotter(formula, accuracy)

    base.start_method()