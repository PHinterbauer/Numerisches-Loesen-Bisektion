import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Base:
    """
    **Base**: A base class for common methods and attributes used.

    Attributes:
        sep_length (int): Length of the separator line.
        methods (list): List of available methods for root finding.
        enable_plot (bool): Indicates if plotting is enabled.
        enable_zoom (bool): Indicates if automatic zooming is enabled for the main plot.
    """
    def __init__(self, 
                formula: str = "", 
                accuracy: float = 0.0, 
                a: float = 0.0, 
                b: float = 0.0, 
                n: float = 0.0, 
                sep_length: int = 100, 
                fig = plt.figure(figsize=(14, 8)), 
                methods: list = ["bisektion", "newton_raphson", "regula_falsi"], 
                enable_plot: bool = False, 
                enable_zoom: bool = False):
        self.formula = formula
        self.accuracy = accuracy
        self.a = a
        self.b = b
        self.n = n
        self.sep_length = sep_length
        self.fig = fig
        self.methods = methods
        self.enable_plot = enable_plot
        self.enable_zoom = enable_zoom

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
        **start_method**: Prompts the user to select a method for root finding and starts the main loop.
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

class Bisektion(Base):
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
        enable_plot (bool): Indicates if plotting is enabled.
        methods (list): List of available methods for root finding.
    """
    # def __init__(self, formula="", accuracy=0.0, fig=plt.figure(figsize=(14, 8)), a=0.0, b=0.0, n=0.0, sep_length=100):
    def __init__(self,
                formula,
                accuracy,
                a, 
                b, 
                n, 
                sep_length, 
                fig, 
                methods, 
                enable_plot, 
                enable_zoom, 
                c: float = 0.0, 
                fa: float = 0.0, 
                fb: float = 0.0, 
                fc: float = 0.0, 
                controll_value: float = 0.0,
                result: dict = {},
                controll: bool = False):
        super().__init__(formula, accuracy, fig, a, b, n, sep_length, fig, methods, enable_plot, enable_zoom)
        self.c = c
        self.fa = fa
        self.fb = fb
        self.fc = fc
        self.controll_value = controll_value
        self.result = result
        self.controll = controll

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
            plotter.initialize_plot(self.a, self.b, self.n)  # Initialize plot
            plotter.update_data(self.a, self.b, self.c, self.fc)  # Update plot data

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
                plotter.update_data(self.a, self.b, self.c, self.fc)  # Update plot data

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
    def __init__(self,
                formula,
                accuracy,
                a, 
                b, 
                n, 
                sep_length, 
                fig, 
                methods, 
                enable_plot, 
                enable_zoom, 
                c, 
                fa, 
                fb, 
                fc, 
                controll_value,
                result,
                controll,
                formula_derivative: str = "",
                derivative: float = 0.0):
        super().__init__(formula, accuracy, a, b, n, sep_length, fig, methods, enable_plot, enable_zoom, c, fa, fb, fc, controll_value, result, controll)
        self.formula_derivative = formula_derivative
        self.derivative = derivative

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
    def __init__(self,
                formula,
                accuracy,
                a, 
                b, 
                n, 
                sep_length, 
                fig, 
                methods, 
                enable_plot, 
                enable_zoom, 
                c, 
                fa, 
                fb, 
                fc, 
                controll_value,
                result,
                controll):
        super().__init__(formula, accuracy, a, b, n, sep_length, fig, methods, enable_plot, enable_zoom, c, fa, fb, fc, controll_value, result, controll)

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation (c) using the Regula Falsi formula.
        """
        # Regula Falsi formula: c = b - (fb * (b - a)) / (fb - fa)
        if self.fb - self.fa == 0:
            raise ValueError("Division by zero in Regula Falsi method.")
        self.c = self.b - (self.fb * (self.b - self.a)) / (self.fb - self.fa)

class Plotter(Base):
    """
    **Plotter**: Handles the visualization of the root-finding process.

    Inherits from:
        Bisektion: Reuses the structure and methods of the bisection method.

    Attributes:
        a_list (list): Stores the 'a' values for each iteration.
        b_list (list): Stores the 'b' values for each iteration.
        c_list (list): Stores the 'c' values (approximations) for each iteration.
        fc_list (list): Stores the function values at 'c' (f(c)) for each iteration.
        c_points_text_str_list (list): Stores the text annotations for 'c' points on the plot.
        enable_zoom (bool): Indicates if automatic zooming is enabled for the main plot.
        fig (matplotlib.figure.Figure): The figure object for the plot.
    """
    def __init__(self, formula="", accuracy=0.0, fig=None, a=0.0, b=0.0, n=0.0, sep_length=100):
        super().__init__(formula, accuracy, fig, a, b, n, sep_length)
        self.a_list: list = []
        self.b_list: list = []
        self.c_list: list = []
        self.fc_list: list = []
        self.c_points_text_str_list: list = []
        self.enable_zoom: bool = False

    def initialize_plot(self, a, b, n):
        """
        **initialize_plot**: Sets up the initial plot for the root-finding process.

        Parameters:
            a (float): Start of the interval.
            b (float): End of the interval.
            n (float): Input number for the formula.
        """
        self.fig.canvas.manager.set_window_title("Nullstellenberechnung durch Iterative-Verfahren")

        # Generate x and y values for the function line
        self.x_func_line = np.linspace(a, b, 500)
        self.y_func_line = [eval(self.formula.format(x=x, n=n)) for x in self.x_func_line]

        # Main plot for root approximation
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

        # Logarithmic plot for |f(c)| per iteration
        self.fc_per_iter_logplot = self.fig.add_subplot(2, 2, 3)
        self.fc_per_iter_logplot_line_list, = self.fc_per_iter_logplot.plot([0], [1e-100], "go-", label="|f(c)|")
        self.fc_per_iter_logplot.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter_logplot.set_title("|f(c)| pro Iteration (Logarithmisch)")
        self.fc_per_iter_logplot.set_xlabel("Iteration")
        self.fc_per_iter_logplot.set_ylabel("|f(c)|")
        self.fc_per_iter_logplot.set_yscale("log")
        self.fc_per_iter_logplot.grid(True)
        self.fc_per_iter_logplot.legend(loc="upper right")

        # Linear plot for f(c) per iteration
        self.fc_per_iter = self.fig.add_subplot(2, 2, 4)
        self.fc_per_iter_line_list, = self.fc_per_iter.plot([], [], "go-", label="f(c)")
        self.fc_per_iter.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter.set_title("f(c) pro Iteration")
        self.fc_per_iter.set_xlabel("Iteration")
        self.fc_per_iter.set_ylabel("f(c)")
        self.fc_per_iter.grid(True)
        self.fc_per_iter.legend(loc="upper right")

    def update_data(self, a, b, c, fc):
        """
        **update_data**: Updates the data lists with new values for each iteration.

        Parameters:
            a (float): Current 'a' value.
            b (float): Current 'b' value.
            c (float): Current 'c' value (approximation).
            fc (float): Current function value at 'c' (f(c)).
        """
        self.a_list.append(a)
        self.b_list.append(b)
        self.c_list.append(c)
        self.fc_list.append(fc)

    def update_plot(self, frame):
        """
        **update_plot**: Updates the plot for each frame during the animation.

        Parameters:
            frame (int): The current frame number in the animation.

        Returns:
            Updated plot elements.
        """
        print("Current a:", self.a_list[frame])
        print("Current b:", self.b_list[frame])
        print("Current c:", self.c_list[frame])
        print("Current fc:", self.fc_list[frame])
        print("Current n:", self.n)
        print("Current iteration:", frame + 1)
        # Clear previous text annotations
        for text_str in self.c_points_text_str_list:
            text_str.remove()
        self.c_points_text_str_list.clear()

        # Update 'a' and 'b' lines aswell as 'c' points
        a_iter = self.a_list[frame]
        b_iter = self.b_list[frame]
        c_point_x_list = self.c_list[:frame + 1]
        c_point_y_list = self.fc_list[:frame + 1]
        self.c_point_list.set_data(c_point_x_list, c_point_y_list)
        self.a_vertical_line.set_xdata([a_iter])
        self.b_vertical_line.set_xdata([b_iter])

        # Adjust zoom if enabled
        if self.enable_zoom:
            approach_to_root_min_width = abs(a_iter - b_iter) # Calculate the width of the interval
            approach_to_root_zoom = approach_to_root_min_width * 0.05 # Add a small margin for zooming
            self.approach_to_root.set_xlim(a_iter - approach_to_root_zoom, b_iter + approach_to_root_zoom)

        # Add text annotations for 'c' points
        c_point_text_y_min, c_point_text_y_max = self.approach_to_root.get_ylim()
        c_point_text_y_offset = (c_point_text_y_max - c_point_text_y_min) * 0.05 # Offset for text placement
        for iter in range(len(c_point_x_list)):
            c_point_text_x_iter = c_point_x_list[iter]
            c_point_text_y_iter = c_point_y_list[iter]
            c_point_text_str_iter = self.approach_to_root.text(
                c_point_text_x_iter, c_point_text_y_iter + c_point_text_y_offset, str(iter + 1), color="red", fontsize=8, ha="center"
            )
            self.c_points_text_str_list.append(c_point_text_str_iter)

        # Update logarithmic plot for |f(c)|
        log_fc_list = np.abs(self.fc_list[:frame + 1]) + 1e-10  # Take absolute values of f(c) up to the current frame and add a small offset to avoid log(0)
        self.fc_per_iter_logplot_line_list.set_data(range(1, frame + 2), log_fc_list) # Use range(1, frame + 2) to match iteration numbers (1-based indexing)
        self.fc_per_iter_logplot.set_xlim(1, len(self.c_list) + 1) # Set x-axis limits to match the number of iterations
        self.fc_per_iter_logplot.set_ylim(min(log_fc_list) * 0.9, max(log_fc_list) * 1.1) # Adjust y-axis limits with a small margin for better visualization

        # Update linear plot for f(c)
        self.fc_per_iter_line_list.set_data(range(1, frame + 2), self.fc_list[:frame + 1]) # Include all f(c) values up to the current frame
        self.fc_per_iter.set_xlim(1, len(self.c_list) + 1) # Set x-axis limits to match the number of iterations
        self.fc_per_iter.set_ylim(min(self.fc_list) * 1.5, max(self.fc_list) * 1.1) # Adjust y-axis limits with a margin to ensure all points are visible

        return self.c_point_list, self.a_vertical_line, self.b_vertical_line, self.fc_per_iter_logplot_line_list, self.fc_per_iter_line_list

    def plot(self):
        """
        **plot**: Creates and displays the animation for the root-finding process.
        """
        fig_anim = animation.FuncAnimation(self.fig, self.update_plot, frames=len(self.c_list), interval=1000, repeat=False)
        plt.tight_layout() # Adjust the layout to prevent overlapping
        plt.show()

if __name__ == "__main__":
    formula = "{x}**2 - {n}" # Formula for root finding (e.g. "np.sqrt({n}) - {x}" / "{x}**2 - self.n")
    accuracy = 0.001  # Desired accuracy (e.g. 1e-50 / 0.001)
    fig = plt.figure(figsize=(14, 8)) # Create a figure for plotting

    base = Base(formula, accuracy, fig)
    bisektion = Bisektion(formula, accuracy, fig)
    newton_raphson = Newton_Raphson(formula, accuracy, fig)
    regula_falsi = Regula_Falsi(formula, accuracy, fig)
    plotter = Plotter(formula, accuracy, fig)

    newton_raphson.formula_derivative = "2 * {x}"  # Derivative of the formula (e.g. "2 * {x}")

    base.start_method()  # Start the method selection process

# TODO:
# update class docstrings 
# move input methods and separator method in extra Class as parent of all classes 
# change call of start_method to use the new class
# nonetype object has no attribute canvas - creates new fig dor every class initialization
# explain comma in initilization of plotter
# rewrite classes and attributes to be inherited correctly