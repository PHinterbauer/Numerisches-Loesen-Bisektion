# ==================================
# Root Finding Methods with Plotting

# Paul Hinterbauer @ TGM Vienna 2025
# ==================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

class Base():
    """
    **Base**: Base class for handling user input and managing settings for numerical methods.

    Attributes:
        formula (str): The mathematical formula to evaluate.
        accuracy (float): The desired accuracy for the calculations.
        a (float): The lower bound of the interval.
        a_initial (float): The initial value of the lower bound.
        b (float): The upper bound of the interval.
        b_initial (float): The initial value of the upper bound.
        n (float): The parameter used in the formula.
        sep_length (int): The length of the separator for console output.
        methods (list): List of available methods for root finding.
        enable_plot (bool): Flag to enable or disable plotting.
        enable_zoom (bool): Flag to enable or disable zooming in the plot.
    """
    def __init__(self, formula: str, accuracy: float, a: float = 0.0, a_initial: float = 0.0, b: float = 0.0, b_initial: float = 0.0, n: float = 0.0, sep_length: int = 100, methods: list = ["bisektion", "newton_raphson", "regula_falsi"], enable_plot: bool = False, enable_zoom: bool = False):
        self.formula = formula
        self.accuracy = accuracy
        self.a = a
        self.a_initial = a_initial
        self.b = b
        self.b_initial = b_initial
        self.n = n
        self.sep_length = sep_length
        self.methods = methods
        self.enable_plot = enable_plot
        self.enable_zoom = enable_zoom

    def separator(self):
        """
        **separator**: Prints a separator line for console output.
        """
        print("-" * self.sep_length)

    def get_plotter_zoom(self):
        """
        **get_plotter_zoom**: Prompts the user to enable or disable zooming in the plot.
        """
        flag = True
        while flag:
            choice = input("Do you want to activate automatic zooming on the main plot? [y/n]: ").lower()
            self.separator()
            if choice in ["y", "yes", "n", "no"]: # Check valid input
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
        """
        **get_plotter**: Prompts the user to enable or disable plotting and handles zoom settings.
        """
        flag = True
        while flag:
            choice = input("Do you want to plot the outcome? [y/n]: ").lower()
            self.separator()
            if choice in ["y", "yes", "n", "no"]: # Check valid input
                if choice in ["y", "yes"]:
                    self.get_plotter_zoom() # Get zoom setting if plotting enabled
                    self.enable_plot = True
                    flag = False
                else:
                    self.enable_plot = False
                    flag = False
            else:
                print("Invalid selection! Please choose [y/n].")
                continue

    def start_method(self):
        """
        **start_method**: Prompts the user to select a numerical method and starts the calculation process.
        """
        self.separator()
        print("Select the method:\n1. Bisection\n2. Newton-Raphson\n3. Regula Falsi")
        flag = True
        while flag:
            try:
                choice = int(input("Enter the number of the method: "))
                self.separator()
                if choice in [1, 2, 3]: # Valid method choices
                    self.get_plotter()
                    # Dynamically call selected method using eval
                    eval(self.methods[choice - 1] + ".main_loop()")
                    if self.enable_plot:
                        plotter.plot()
                    flag = False
                else:
                    print("Invalid selection! Please choose 1, 2, or 3.")
                    continue
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_n(self):
        """
        **get_n**: Prompts the user to input the parameter 'n' for the formula.
        """
        if "{n}" in self.formula: # Check if formula uses n parameter
            flag = True
            while flag:
                try:
                    self.n = float(input(f'Enter n for: '))
                    flag = False
                except ValueError:
                    print("Please enter a valid number!")
                    continue
        else:
            self.formula += " + {n}" # Add n parameter if not present using default value of 0 so it doesn't affect the formula

    def get_interval_a(self):
        """
        **get_interval_a**: Prompts the user to input the lower bound of the interval.
        """
        flag = True
        while flag:
            try:
                self.a = self.a_initial = float(input(f'Enter a for: '))
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

    def get_interval_b(self):
        """
        **get_interval_b**: Prompts the user to input the upper bound of the interval.
        """
        flag = True
        while flag:
            try:
                self.b = self.b_initial = float(input(f'Enter b for: '))
                self.separator()
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

class BaseCalculations(Base):
    """
    **BaseCalculations**: Base class for performing calculations related to root-finding methods.

    Inherits from:
        Base: Provides user input handling and settings management.

    Attributes:
        base (Base): The base class instance.
        c (float): The midpoint or calculated root.
        fa (float): The function value at 'a'.
        fb (float): The function value at 'b'.
        fc (float): The function value at 'c'.
        control_value (float): The product of 'fa' and 'fb'.
        result (dict): Stores the result of the calculation.
        control (bool): Indicates if the interval contains a root.
    """
    def __init__(self, base: Base, iteration_speed: float, c: float = 0.0, fa: float = 0.0, fb: float = 0.0, fc: float = 0.0, control_value: float = 0.0, result: dict = {}, control: bool = False):
        self.base = base
        self.iteration_speed = iteration_speed
        self.c = c
        self.fa = fa
        self.fb = fb
        self.fc = fc
        self.control_value = control_value
        self.result = result
        self.control = control

    # Delegate attribute access to base class
    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def calculate_fa(self):
        """
        **calculate_fa**: Calculates the function value at 'a'.
        """
        self.fa = eval(self.formula.format(x=self.a, n=self.n))

    def calculate_fb(self):
        """
        **calculate_fb**: Calculates the function value at 'b'.
        """
        self.fb = eval(self.formula.format(x=self.b, n=self.n))

    def calculate_c(self):
        """
        **calculate_c**: Abstract method to calculate the value of 'c'. Must be implemented in subclasses.
        """
        raise NotImplementedError("Die Methode 'calculate_c' muss in der Unterklasse überschrieben werden.")

    def calculate_fc(self):
        """
        **calculate_fc**: Calculates the function value at 'c'.
        """
        self.fc = eval(self.formula.format(x=self.c, n=self.n))

    def check_control(self):
        """
        **check_control**: Checks if the interval contains a root by evaluating the product of 'fa' and 'fb'.
        """
        self.control_value = self.fa * self.fb
        self.control = self.control_value < 0 # Root exists if signs differ

    def switch_interval(self):
        """
        **switch_interval**: Updates the interval based on the sign of the function values.
        """
        if self.fa * self.fc < 0: # Root in left subinterval
            self.b = self.c
        else: # Root in right subinterval
            self.a = self.c 

    def calc_separator(self, text_1: str, header: str, text_2: str = "", text_3: str = ""):
        """
        **calc_separator**: Calculates the padding for formatted console output.

        Parameters:
            text_1 (str): The first line of text.
            header (str): The header text.
            text_2 (str): The second line of text (optional).
            text_3 (str): The third line of text (optional).

        Returns:
            left_side_length (int): Left padding.
            right_side_length (int): Right padding.
            total_length (int): Total length of the separator.
        """
        total_length = max(len(text_1), len(text_2), len(text_3))
        side_length = (total_length - len(header)) // 2 # Center header

        left_side_length = side_length
        right_side_length = total_length - len(header) - left_side_length

        return left_side_length, right_side_length, total_length

    def main_loop(self):
        """
        **main_loop**: Executes the main iterative process for root-finding.
        """ 
        # Initialize variables
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
            plotter.update_data(self.a, self.b, self.c, self.fc, self.n, self.a_initial, self.b_initial)

        # Format and print iteration info
        left_side_length, right_side_length, total_length = self.calc_separator(
            text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
            text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
            text_3=f'control: {self.control}, control value: {self.control_value}',
            header=f' Iteration: {iteration} ')
        print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
        print(f'a: {self.a}, b: {self.b}, c: {self.c}')
        print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
        print(f'control: {self.control}, control value: {self.control_value}')
        print("=" * total_length)

        self.switch_interval() # Switch interval for next iteration

        time.sleep(self.iteration_speed) # Pause for visualization

        # Continue until fc reaches desired accuracy or root found
        while np.floor(abs(self.fc) * 10 ** (abs(int(np.log10(self.accuracy))) - 1)) != 0 and self.control: # Using floor to avoid floating point errors and * 10 ** (abs(int(np.log10(self.accuracy))) - 1) to get the correct number of decimal places
            # Calculate next iteration
            iteration += 1
            self.calculate_fa()
            self.calculate_fb()
            self.calculate_c()
            self.calculate_fc()
            self.check_control()

            if self.enable_plot:
                plotter.update_data(self.a, self.b, self.c, self.fc, self.n, self.a_initial, self.b_initial)

            # Format and print iteration info
            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                text_2=f'a: {self.a}, b: {self.b}, c: {self.c}',
                text_3=f'control: {self.control}, control value: {self.control_value}',
                header=f' Iteration: {iteration} ')
            print("=" * left_side_length, "Iteration: {iteration}".format(iteration=iteration), "=" * right_side_length)
            print(f'a: {self.a}, b: {self.b}, c: {self.c}')
            print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
            print(f'control: {self.control}, control value: {self.control_value}')
            print("=" * total_length)

            self.switch_interval() # Switch interval for next iteration

            time.sleep(self.iteration_speed) # Pause for visualization

            # Check for exact root found
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
            self.result["c"] = self.c # Store final result

        # Print final results
        if self.result:
            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                header=f' Result ')
            print("=" * left_side_length, "Result", "=" * right_side_length)
            for key, value in self.result.items():
                result_text = f'{key}: {value}'
                left_padding, right_padding, _ = self.calc_separator(text_1="=" * total_length, header=result_text)
                print(" " * left_padding + result_text + " " * right_padding)
            print("=" * total_length)
        else:
            print("No solution found!")

class Bisektion(BaseCalculations):
    """
    **Bisektion**: Implements the bisection method for root-finding.

    Inherits from:
        BaseCalculations: Provides the base functionality for calculations.

    Methods:
        calculate_c: Calculates the midpoint of the interval.
    """
    def __init__(self, base: BaseCalculations):
        self.base = base
    
    # Delegate attribute access to base class
    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def calculate_c(self):
        """
        **calculate_c**: Calculates the midpoint of the interval.
        """
        self.c = (self.a + self.b) / 2 

class NewtonRaphson(BaseCalculations):
    """
    **NewtonRaphson**: Implements the Newton-Raphson method for root-finding.

    Inherits from:
        BaseCalculations: Provides the base functionality for calculations.

    Attributes:
        formula_derivative (str): The derivative of the formula.
        c_derivative (float): The derivative value at 'c'.
    """
    def __init__(self, base: BaseCalculations, formula_derivative: str, c_derivative: float = 0.0):
        self.base = base
        self.formula_derivative = formula_derivative
        self.c_derivative = c_derivative

    # Delegate attribute access to base class
    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation of the root using the Newton-Raphson formula.
        """
        self.derivative = eval(self.formula_derivative.format(x=self.c)) # Calculate derivative
        try:
            self.c = self.c - self.fc / self.derivative # Newton-Raphson formula
        except ZeroDivisionError:
            self.separator()
            print("Division by zero! Offsetting by 1e-10")
            self.separator()
            self.c = self.c - self.fc / (self.derivative + 1e-10) # Handle division by zero

class RegulaFalsi(BaseCalculations):
    """
    **RegulaFalsi**: Implements the Regula Falsi method for root-finding.

    Inherits from:
        BaseCalculations: Provides the base functionality for calculations.

    Methods:
        calculate_c: Calculates the next approximation of the root using the Regula Falsi formula.
    """
    def __init__(self, base: BaseCalculations):
        self.base = base

    # Delegate attribute access to base class
    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation of the root using the Regula Falsi formula.
        """
        self.c = (self.a * self.fb - self.b * self.fa) / (self.fb - self.fa) # Regula Falsi formula

class Plotter(BaseCalculations):
    """
    **Plotter**: Handles the visualization of the root-finding process.

    Inherits from:
        BaseCalculations: Provides the base functionality for calculations.

    Attributes:
        plot_speed (int): Speed of the plot animation.
        a_list (list): List of 'a' values over iterations.
        b_list (list): List of 'b' values over iterations.
        c_list (list): List of 'c' values over iterations.
        fc_list (list): List of 'fc' values over iterations.
        c_points_text_str_list (list): List of text annotations for 'c' points.
    """
    def __init__(self, base: BaseCalculations, plot_speed: int, a_list: list = [], b_list: list = [], c_list: list = [], fc_list: list = [], c_points_text_str_list: list = []):
        self.base = base
        self.plot_speed = plot_speed
        self.a_list = a_list
        self.b_list = b_list
        self.c_list = c_list
        self.fc_list = fc_list
        self.c_points_text_str_list = c_points_text_str_list

    # Delegate attribute access to base class
    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def initialize_plot(self):
        """
        **initialize_plot**: Sets up the plot layout and initializes plot elements.
        """
        # Create figure and set title
        self.fig = plt.figure(figsize=(14, 8))
        self.fig.canvas.manager.set_window_title("Nullstellenberechnung durch Iterative-Verfahren")

        # Create x values and calculate y values for function plot
        self.x_func_line = np.linspace(self.a_initial, self.b_initial, 500)
        self.y_func_line = [eval(self.formula.format(x=x, n=self.n)) for x in self.x_func_line]

        # Main subplot for function plot
        self.approach_to_root = self.fig.add_subplot(2, 1, 1)
        self.func_line_list, = self.approach_to_root.plot(self.x_func_line, self.y_func_line, label=self.formula.replace("{", "").replace("}", ""), color="blue") # Unpacks single-element tuple to get the line object directly
        self.c_point_list, = self.approach_to_root.plot([], [], "ro", label="Geschätzte Nullstelle") # Unpacks single-element tuple to get the red circle (ro) object directly
        self.a_vertical_line = self.approach_to_root.axvline(0, color="green", linestyle="--", label="a")
        self.b_vertical_line = self.approach_to_root.axvline(0, color="purple", linestyle="--", label="b")
        self.approach_to_root.axhline(0, color="gray", linestyle="--")
        self.approach_to_root.set_title("Annäherung an die Nullstelle")
        self.approach_to_root.set_xlabel("c")
        self.approach_to_root.set_ylabel("f(c)")
        self.approach_to_root.grid(True)
        self.approach_to_root.legend(loc="upper left")

        # Log plot of |f(c)| over iterations 
        self.fc_per_iter_logplot = self.fig.add_subplot(2, 2, 3)
        self.fc_per_iter_logplot_line_list, = self.fc_per_iter_logplot.plot([0], [1e-100], "go-", label="|f(c)|") # Unpacks single-element tuple to get the line object directly
        self.fc_per_iter_logplot.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter_logplot.set_title("|f(c)| pro Iteration (Logarithmisch)")
        self.fc_per_iter_logplot.set_xlabel("Iteration")
        self.fc_per_iter_logplot.set_ylabel("|f(c)|")
        self.fc_per_iter_logplot.set_yscale("log") # Logarithmic y-axis
        self.fc_per_iter_logplot.grid(True)
        self.fc_per_iter_logplot.legend(loc="upper right")

        # Linear plot of f(c) over iterations
        self.fc_per_iter = self.fig.add_subplot(2, 2, 4)
        self.fc_per_iter_line_list, = self.fc_per_iter.plot([], [], "go-", label="f(c)") # Unpacks single-element tuple to get the line object directly
        self.fc_per_iter.axhline(0, color="gray", linestyle="--")
        self.fc_per_iter.set_title("f(c) pro Iteration")
        self.fc_per_iter.set_xlabel("Iteration")
        self.fc_per_iter.set_ylabel("f(c)")
        self.fc_per_iter.grid(True)
        self.fc_per_iter.legend(loc="upper right")

    def update_data(self, a ,b, c, fc, n, a_initial, b_initial):
        """
        **update_data**: Updates the data for plotting.

        Parameters:
            a (float): The current value of 'a'.
            b (float): The current value of 'b'.
            c (float): The current value of 'c'.
            fc (float): The function value at 'c'.
            n (float): The parameter used in the formula.
            a_initial (float): The initial value of 'a'.
            b_initial (float): The initial value of 'b'.
        """
        self.n = n
        self.a_initial = a_initial
        self.b_initial = b_initial
        self.a_list.append(a)
        self.b_list.append(b)
        self.c_list.append(c)
        self.fc_list.append(fc)

    def update_plot(self, frame):
        """
        **update_plot**: Updates the plot elements for each frame in the animation.

        Parameters:
            frame (int): The current frame index.

        Returns:
            c_point_list (list): Updated data for 'c' points.
            a_vertical_line (list): Updated vertical line for 'a'.
            b_vertical_line (list): Updated vertical line for 'b'.
            fc_per_iter_logplot_line_list (list): Updated log plot data for |f(c)|.
            fc_per_iter_line_list (list): Updated plot data for f(c).
        """
        # Clear previous text annotations
        for text_str in self.c_points_text_str_list:
            text_str.remove()
        self.c_points_text_str_list.clear()

        # Get data for current frame
        a_iter = self.a_list[frame]
        b_iter = self.b_list[frame]
        c_point_x_list = self.c_list[:frame + 1]
        c_point_y_list = self.fc_list[:frame + 1]

        # Update plot elements with current data
        self.c_point_list.set_data(c_point_x_list, c_point_y_list)
        self.a_vertical_line.set_xdata([a_iter])
        self.b_vertical_line.set_xdata([b_iter])

        # Auto-zoom if enabled
        if self.enable_zoom:
            approach_to_root_min_width = abs(a_iter - b_iter)
            approach_to_root_zoom = approach_to_root_min_width * 0.05
            self.approach_to_root.set_xlim(a_iter - approach_to_root_zoom, b_iter + approach_to_root_zoom)

        # Add iteration number annotations
        c_point_text_y_min, c_point_text_y_max = self.approach_to_root.get_ylim()
        c_point_text_y_offset = (c_point_text_y_max - c_point_text_y_min) * 0.05
        for iter in range(len(c_point_x_list)):
            c_point_text_x_iter = c_point_x_list[iter]
            c_point_text_y_iter = c_point_y_list[iter]
            c_point_text_str_iter = self.approach_to_root.text(c_point_text_x_iter, c_point_text_y_iter + c_point_text_y_offset, str(iter + 1), color="red", fontsize=8, ha="center")
            self.c_points_text_str_list.append(c_point_text_str_iter)

        # Update log plot
        log_fc_list = np.abs(self.fc_list[:frame + 1]) + 1e-10 # Avoid log(0)
        self.fc_per_iter_logplot_line_list.set_data(range(1, frame + 2), log_fc_list)
        self.fc_per_iter_logplot.set_xlim(1, len(self.c_list) + 1)
        self.fc_per_iter_logplot.set_ylim(min(log_fc_list) * 0.9, max(log_fc_list) * 1.1)

        # Update linear plot
        self.fc_per_iter_line_list.set_data(range(1, frame + 2), self.fc_list[:frame + 1])
        self.fc_per_iter.set_xlim(1, len(self.c_list) + 1)
        self.fc_per_iter.set_ylim(min(self.fc_list) * 1.5, max(self.fc_list) * 1.1)

        return self.c_point_list, self.a_vertical_line, self.b_vertical_line, self.fc_per_iter_logplot_line_list, self.fc_per_iter_line_list

    def plot(self):
        """
        **plot**: Initializes and displays the plot animation.
        """
        self.initialize_plot()
        # Create animation with update_plot function
        anim1 = animation.FuncAnimation(self.fig, self.update_plot, frames=len(self.c_list), interval=self.plot_speed, repeat=False)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    # Example formula and its derivative aswell as accuracy and speed settings
    formula = "{x}**2 - {n} " # e.g. {x}**2 - {n} / 2*{x} + {x}**2 + 3*{x}**3 - {x}**4
    formula_derivative = "2 * {x}" # e.g. 2 * {x} / 2 + 2*{x} + 9*{x}**2 - 4*{x}**3
    accuracy = 0.001 # e.g. 1e-50
    plot_speed = 1000 # ms per frame
    iteration_speed = 0.2 # s per iteration

    # Initialize all components
    base = Base(formula, accuracy)
    base_calculations = BaseCalculations(base, iteration_speed)
    bisektion = Bisektion(base_calculations)
    newton_raphson = NewtonRaphson(base_calculations, formula_derivative)
    regula_falsi = RegulaFalsi(base_calculations)
    plotter = Plotter(base_calculations, plot_speed)

    # Start the selected method
    base.start_method()
