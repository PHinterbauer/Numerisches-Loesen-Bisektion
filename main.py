import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Bisektion():
    """
    **Bisektion**: Implements the bisection method for root finding.

    Attributes:
        a (float): Start of interval.
        b (float): End of interval.
        c (float): Midpoint of interval.
        n (float): Input number for the formula.
        fa (float): Function value at a.
        fb (float): Function value at b.
        fc (float): Function value at c.
        accuracy (float): Desired accuracy for the result.
        controll_value (float): Product of fa and fb.
        result (dict): Stores the final result.
        formula (str): Mathematical formula as a string.
        controll (bool): Indicates if the interval is valid.
        methods (list): List of available methods for root finding.
    """
    def __init__(self, formula: str = "", accuracy: float = 0.0):
        self.a: float = 0.0
        self.b: float = 0.0
        self.c: float = 0.0
        self.n: float = 0.0
        self.fa: float = 0.0
        self.fb: float = 0.0
        self.fc: float = 0.0
        self.accuracy = accuracy
        self.controll_value: float = 0.0
        self.sep_length: int = 50
        self.result: dict = {}
        self.formula = formula
        self.controll: bool = False
        self.enable_plot: bool = False
        self.enable_zoom: bool = False
        self.methods: list = ["bisektion", "newton_raphson", "regula_falsi"]

    def separator(self):
        """
        **separator**: Prints a separator line of specified length.
        """
        print("-" * self.sep_length)

    def get_plotter_zoom(self):
        """
        **get_plotter_zoom**: Prompts the user if they want automatic zooming on the main plot.
        """
        flag = True
        while flag:
            choice = input("Do you want to activate automatic zooming on the main plot? [y/n]: ").lower() # Input for zooming
            self.separator()  # Print separator
            if choice in ["y", "yes", "n", "no"]:
                if choice in ["y", "yes"]:
                    self.enable_zoom = True # Set zoom flag to True
                    flag = False
                else:
                    self.enable_zoom = False # Set zoom flag to False
                    flag = False
            else:
                print("Invalid selection! Please choose [y/n].")
                continue

    def get_plotter(self):
        """
        **get_plotter**: Prompts the user if they want to plot the outcome.
        """
        flag = True
        while flag:
            choice = input("Do you want to plot the outcome? [y/n]: ").lower()  # Input for plotting
            self.separator()  # Print separator
            if choice in ["y", "yes", "n", "no"]:
                if choice in ["y", "yes"]:
                    self.get_plotter_zoom()
                    self.enable_plot = True # Set plot flag to True
                    flag = False
                else:
                    self.enable_plot = False # Set plot flag to False
                    flag = False
            else:
                print("Invalid selection! Please choose [y/n].")
                continue

    def start_method(self):
        """
        **get_method**: Prompts the user to select a method for root finding.
        """
        self.separator()  # Print separator
        print("Select the method:\n1. Bisection\n2. Newton-Raphson\n3. Regula Falsi")
        flag = True
        while flag:
            try:
                choice = int(input("Enter the number of the method: "))  # Input for method selection
                self.separator()  # Print separator
                if choice in [1, 2, 3]:
                    self.get_plotter()
                    eval(self.methods[choice - 1] + ".main_loop()")  # Call the selected method's main loop
                    if self.enable_plot:
                        plotter.plot()  # Call the plot method if selected
                    flag = False
                else:
                    print("Invalid selection! Please choose 1, 2, or 3.")
                    continue
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_n(self):
        """
        **get_n**: Prompts the user to input the value of n.
        """
        flag = True
        while flag:
            try:
                self.n = float(input(f'Enter n for: '))  # Input for n
                flag = False
            except ValueError:
                print("Please enter a valid number!")  # Handle invalid input
                continue

    def get_interval_a(self):
        """
        **get_interval_a**: Prompts the user to input the start of the interval (a).
        """
        flag = True
        while flag:
            try:
                self.a = float(input(f'Enter a for: '))  # Input for interval start
                flag = False
            except ValueError:
                print("Please enter a valid number!")  # Handle invalid input
                continue

    def get_interval_b(self):
        """
        **get_interval_b**: Prompts the user to input the end of the interval (b).
        """
        flag = True
        while flag:
            try:
                self.b = float(input(f'Enter b for: '))  # Input for interval end
                self.separator()  # Print separator
                flag = False
            except ValueError:
                print("Please enter a valid number!")  # Handle invalid input
                continue

    def calculate_fa(self):
        """
        **calculate_fa**: Calculates the function value at a (fa).
        """
        self.fa = eval(self.formula.format(x=self.a, n=self.n))  # Evaluate formula at a

    def calculate_fb(self):
        """
        **calculate_fb**: Calculates the function value at b (fb).
        """
        self.fb = eval(self.formula.format(x=self.b, n=self.n))  # Evaluate formula at b

    def calculate_c(self):
        """
        **calculate_c**: Calculates the midpoint of the interval (c).
        """
        self.c = (self.a + self.b) / 2  # Midpoint formula

    def calculate_fc(self):
        """
        **calculate_fc**: Calculates the function value at c (fc).
        """
        self.fc = eval(self.formula.format(x=self.c, n=self.n))  # Evaluate formula at c

    def check_controll(self):
        """
        **check_controll**: Checks if the interval is valid by evaluating fa * fb.
        """
        self.controll_value = self.fa * self.fb  # Product of fa and fb
        self.controll = self.controll_value < 0 # Valid if product is negative

    def switch_interval(self):
        """
        **switch_interval**: Updates the interval based on the sign of fa * fc.
        """
        if self.fa * self.fc < 0:
            self.b = self.c  # Update b if fa * fc < 0
        else:
            self.a = self.c  # Update a otherwise

    def calc_separator(self, text_1: str, header: str, text_2: str = "", text_3: str = ""):
        """
        **calc_separator**: Calculates padding for formatted output.

        Parameters:
            text_1 (str): The first text to align (usually the longest).
            header (str): The header text.
            text_2 (str, optional): A second text to consider for length.
            text_3 (str, optional): A third text to consider for length.

        Returns:
            left_side_length (int), right_side_length (int), and total length (int).
        """
        # Determine the longest text among the provided options
        total_length = max(len(text_1), len(text_2), len(text_3))
        side_length = (total_length - len(header)) // 2  # Calculate side padding

        left_side_length = side_length
        right_side_length = total_length - len(header) - left_side_length

        return left_side_length, right_side_length, total_length

    def main_loop(self):
        """
        **main_loop**: Executes the bisection method until the desired accuracy is achieved.
        """
        iteration = 1
        self.result = {}
        self.get_n()  # Get input for n
        self.get_interval_a()  # Get interval start
        self.get_interval_b()  # Get interval end
        self.calculate_fa()  # Calculate fa
        self.calculate_fb()  # Calculate fb
        self.calculate_c()  # Calculate midpoint c
        self.calculate_fc()  # Calculate fc
        self.check_controll()  # Check if interval is valid

        if self.enable_plot:
            plotter.update_data()  # Update plot data

        # Print initial iteration
        left_side_length, right_side_length, total_length = self.calc_separator(
            text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
            text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
            text_3=f'controll: {self.controll}, controll value: {self.controll_value}',
            header=f' Iteration: {iteration} '
        )
        print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
        print(f'a: {self.a}, b: {self.b}, c: {self.c}')
        print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
        print(f'controll: {self.controll}, controll value: {self.controll_value}')
        print("=" * total_length)

        self.switch_interval()  # Update interval

        # Perform iterations
        while np.floor(abs(self.fc) * 10 ** (abs(int(np.log10(self.accuracy))) - 1)) != 0 and self.controll:
            iteration += 1
            self.calculate_fa()  # Update fa
            self.calculate_fb()  # Update fb
            self.calculate_c()  # Update c
            self.calculate_fc()  # Update fc
            self.check_controll()  # Recheck interval validity

            if self.enable_plot:
                plotter.update_data()  # Update plot data

            # Print iterations
            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
                text_3=f'controll: {self.controll}, controll value: {self.controll_value}',
                header=f' Iteration: {iteration} '
            )
            print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
            print(f'a: {self.a}, b: {self.b}, c: {self.c}')
            print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
            print(f'controll: {self.controll}, controll value: {self.controll_value}')
            print("=" * total_length)

            self.switch_interval()  # Update interval

            # Check for exact solutions
            if self.fc == 0:
                self.c = self.c  # Exact solution at c
                break
            elif self.fa == 0:
                self.c = self.a  # Exact solution at a
                break
            elif self.fb == 0:
                self.c = self.b  # Exact solution at b
                break

        # Check if the loop ended due to a valid solution or not
        if not self.controll:
            print("No solution found!") # No valid solution
        else:
            self.result["c"] = self.c

        # Print the result
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
            print("No solution found!")  # No valid solution

class Newton_Raphson(Bisektion):
    """
    **Newton_Raphson**: Implements the Newton-Raphson method for root finding.

    Inherits from:
        Bisektion: Reuses the structure and methods of the bisection method.

    Attributes:
        formula (str): Mathematical formula as a string.
        formula_derivative (str): The derivative of the mathematical formula as a string.
        derivative (float): The derivative value at the current approximation (c).
    """
    def __init__(self, formula: str = "", accuracy: float = 0.0):
        super().__init__(formula=formula, accuracy=accuracy)
        self.formula_derivative: str = ""
        self.derivative: float = 0.0

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation (c) using the Newton-Raphson formula.
        """
        # Ensure c is not initialized to a value that causes the derivative to be zero
        if self.c == 0:
            if self.a == 0:
                self.c = self.b
            else:
                self.c = self.a

        # Newton-Raphson formula: c = c - f(c) / f'(c)
        self.derivative = eval(self.formula_derivative.format(x=self.c))  # Evaluate derivative at c
        if self.derivative == 0:
            raise ValueError("Derivative is zero. Newton-Raphson method fails.")
        self.c = self.c - self.fc / self.derivative  # Update c

class Regula_Falsi(Bisektion):
    """
    **Regula_Falsi**: Implements the Regula Falsi (False Position) method for root finding.

    Inherits from:
        Bisektion: Reuses the structure and methods of the bisection method.
    """
    def __init__(self, formula: str = "", accuracy: float = 0.0):
        super().__init__(formula=formula, accuracy=accuracy)

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation (c) using the Regula Falsi formula.
        """
        # Regula Falsi formula: c = b - (fb * (b - a)) / (fb - fa)
        if self.fb - self.fa == 0:
            raise ValueError("Division by zero in Regula Falsi method.")
        self.c = self.b - (self.fb * (self.b - self.a)) / (self.fb - self.fa)

class Plotter(Bisektion):
    def __init__(self):
        super().__init__()
        self.a_list: list = []
        self.b_list: list = []
        self.c_list: list = []
        self.fc_list: list = []
        self.c_points_text_str_list: list = []

    def initialize_plot(self):
        self.fig = plt.figure(figsize=(14, 8))
        self.fig.canvas.manager.set_window_title("Nullstellenberechnung durch Iterative-Verfahren")  # Set the figure's window title
        self.x_func_line = np.linspace(self.a, self.b, 500) # Generate x values for plotting the function line
        self.y_func_line = [eval(self.formula.format(x=x, n=self.n)) for x in self.x_func_line]  # Generate y values using the formula for plotting the function line

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
        self.fc_per_iter_logplot_line_list, = self.fc_per_iter_logplot.plot([], [], "go-", label="|f(c)|") # Initialize with a small value to always get positive values for log
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
        
    def update_data(self):
        self.a_list.append(self.a)
        self.b_list.append(self.b)
        self.c_list.append(self.c)
        self.fc_list.append(self.fc)

    def update_plot(self, frame):
        for text_str in self.c_points_text_str_list:
            text_str.remove()
        self.c_points_text_str_list.clear()

        a_iter = self.a_list[frame]
        b_iter = self.b_list[frame]
        c_iter = self.c_list[frame]
        fc_iter = self.fc_list[frame]
        self.c_point_list.set_data([c_iter], [fc_iter]) # possible error
        self.a_vertical_line.set_xdata([a_iter])
        self.b_vertical_line.set_xdata([b_iter])

        if self.enable_zoom:
            approach_to_root_min_width = a_iter - b_iter
            approach_to_root_zoom = approach_to_root_min_width * 0.05
            self.approach_to_root.set_xlim(a_iter - approach_to_root_zoom, b_iter + approach_to_root_zoom)

        for iter in range(frame + 1):
            c_point_text_x_iter = self.c_list[iter]
            c_point_text_y_iter = self.fc_list[iter]
            c_point_text_str_iter = self.approach_to_root.text(c_point_text_x_iter, c_point_text_y_iter + 20, str(iter), color="red", fontsize=8, ha='center')
            self.c_points_text_str_list.append(c_point_text_str_iter)

        self.fc_per_iter_logplot_line_list.set_data(range(frame + 1), self.fc_list[:frame + 1])
        self.fc_per_iter_logplot.set_xlim(0, frame + 1)
        safe_fc_vals = np.abs(self.fc_list) + 1e-10  # Avoid log(0) by adding a small value
        self.fc_per_iter_logplot.set_ylim(min(safe_fc_vals) * 0.9, max(safe_fc_vals) * 1.1)

        self.fc_per_iter_line_list.set_data(range(frame + 1), self.fc_list[:frame + 1])
        self.fc_per_iter.set_xlim(0, frame + 1)
        ymin = min(self.fc_list)
        ymax = max(self.fc_list)
        if abs(ymax - ymin) < 1e-8:
            ymin -= 1
            ymax += 1
        self.fc_per_iter.set_ylim(ymin, ymax)
        self.fc_per_iter.set_ylim(ymin * 1.1, ymax * 1.1)

        return self.c_point_list, self.a_vertical_line, self.b_vertical_line, self.fc_per_iter_logplot_line_list, self.fc_per_iter_line_list

    def plot(self):
        self.initialize_plot()
        fig_anim = animation.FuncAnimation(self.fig, self.update_plot, frames=len(self.c_list), interval=1000, repeat=False)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    bisektion = Bisektion()
    newton_raphson = Newton_Raphson()
    regula_falsi = Regula_Falsi()
    plotter = Plotter()

    bisektion.formula = "{x}**2 - {n}" # Formula for root finding (e.g. "np.sqrt({n}) - {x}" / "{x}**2 - self.n")
    plotter.formula = bisektion.formula # Set formula for plotter

    newton_raphson.formula_derivative = "2 * {x}"  # Derivative of the formula (e.g. "2 * {x}")

    bisektion.accuracy = 0.001  # Desired accuracy (e.g. 1e-50 / 0.001)

    bisektion.start_method()  # Start the method selection process

# write plotter class from TEST_16.py
# fix no positive values for log warning
# add comments and docstrings to every function that is missing them
# -maybe- write gui class for options