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
    def __init__(self):
        self.a: float = 0.0
        self.b: float = 0.0
        self.c: float = 0.0
        self.n: float = 0.0
        self.fa: float = 0.0
        self.fb: float = 0.0
        self.fc: float = 0.0
        self.accuracy: float = 0.0
        self.controll_value: float = 0.0
        self.sep_length: int = 50
        self.result: dict = {}
        self.formula: str = ""
        self.controll: bool = False
        self.enable_plot: bool = False
        self.methods: list = ["bisektion", "newton_raphson", "regula_falsi"]

    def separator(self):
        """
        **separator**: Prints a separator line of specified length.
        """
        print("-" * self.sep_length)

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
                    self.enable_plot = True # Set plot flag to True
                    flag = False
                else:
                    self.enable_plot = False
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
            plotter.update_data(iteration, self.fc, self.c, (self.a, self.b))  # Update plot data

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
                plotter.update_data(iteration, self.fc, self.c, (self.a, self.b))  # Update plot data


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
    def __init__(self):
        super().__init__()
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
    def __init__(self):
        super().__init__()

    def calculate_c(self):
        """
        **calculate_c**: Calculates the next approximation (c) using the Regula Falsi formula.
        """
        # Regula Falsi formula: c = b - (fb * (b - a)) / (fb - fa)
        if self.fb - self.fa == 0:
            raise ValueError("Division by zero in Regula Falsi method.")
        self.c = self.b - (self.fb * (self.b - self.a)) / (self.fb - self.fa)

# class Plotter(Bisektion):
#     def __init__(self):
#         super().__init__()
#         self.iterations = []
#         self.accuracies = []
#         self.solutions = []
#         self.intervals = []

#     def update_data(self, iteration, accuracy, solution, interval):
#         """
#         Updates the data for plotting.
#         """
#         self.iterations.append(iteration)
#         self.accuracies.append(accuracy)
#         self.solutions.append(solution)
#         self.intervals.append(interval)

#     def plot(self):
#         """
#         Creates animated plots for Aufgabe 7.
#         """
#         fig, axes = plt.subplots(3, 1, figsize=(12, 12))  # Create 3 subplots
#         ax1, ax2, ax3 = axes

#         # Subplot 1: Interval shrinking visualization
#         ax1.set_title("Interval Shrinking Over Iterations")
#         ax1.set_xlabel("Iteration")
#         ax1.set_ylabel("Interval Size")
#         line1, = ax1.plot([], [], 'b-', label="Interval Size")
#         ax1.legend()
#         ax1.grid()

#         # Subplot 2: Accuracy convergence
#         ax2.set_title("Accuracy Convergence (|f(c)|)")
#         ax2.set_xlabel("Iteration")
#         ax2.set_ylabel("Accuracy (|f(c)|)")
#         line2, = ax2.plot([], [], 'g-', label="Accuracy")
#         ax2.legend()
#         ax2.grid()

#         # Subplot 3: Solution convergence
#         ax3.set_title("Solution Convergence (c)")
#         ax3.set_xlabel("Iteration")
#         ax3.set_ylabel("Solution (c)")
#         line3, = ax3.plot([], [], 'r-', label="Solution")
#         ax3.legend()
#         ax3.grid()

#         # Data for animations
#         interval_sizes = [abs(b - a) for a, b in self.intervals]

#         def init():
#             line1.set_data([], [])
#             line2.set_data([], [])
#             line3.set_data([], [])
#             return line1, line2, line3

#         def update(frame):
#             x_data = self.iterations[:frame + 1]
#             y1_data = interval_sizes[:frame + 1]
#             y2_data = self.accuracies[:frame + 1]
#             y3_data = self.solutions[:frame + 1]

#             line1.set_data(x_data, y1_data)
#             line2.set_data(x_data, y2_data)
#             line3.set_data(x_data, y3_data)

#             ax1.relim()
#             ax1.autoscale_view()
#             ax2.relim()
#             ax2.autoscale_view()
#             ax3.relim()
#             ax3.autoscale_view()

#             return line1, line2, line3

#         # Create the animation with a slower interval (e.g., 500ms per frame)
#         ani = animation.FuncAnimation(fig, update, frames=len(self.iterations), init_func=init, interval=500, repeat=False)

#         plt.tight_layout()
#         plt.show()

class Plotter(Bisektion):
    """
    **Plotter**: Handles data collection and visualization for the root-finding process.

    Inherits from:
        Bisektion: Reuses the structure and methods of the bisection method.

    Attributes:
        iterations (list): Stores the iteration numbers.
        accuracies (list): Stores the accuracy values (f(c)) for each iteration.
        solutions (list): Stores the solution values (c) for each iteration.
        intervals (list): Stores the interval bounds (a, b) for each iteration.
    """
    def __init__(self):
        super().__init__()
        self.iterations = []
        self.accuracies = []
        self.solutions = []
        self.intervals = []

    def update_data(self, iteration, accuracy, solution, interval):
        """
        **update_data**: Updates the data for plotting after each iteration.

        Parameters:
            iteration (int): The current iteration number.
            accuracy (float): The accuracy value (f(c)) for the current iteration.
            solution (float): The solution value (c) for the current iteration.
            interval (tuple): The interval bounds (a, b) for the current iteration.
        """
        self.iterations.append(iteration)  # Append the iteration number
        self.accuracies.append(accuracy)  # Append the accuracy value
        self.solutions.append(solution)  # Append the solution value
        self.intervals.append(interval)  # Append the interval bounds

if __name__ == "__main__":
    bisektion = Bisektion()
    newton_raphson = Newton_Raphson()
    regula_falsi = Regula_Falsi()
    plotter = Plotter()

    bisektion.formula = "{x}**2 - {n}" # Formula for root finding (e.g. "np.sqrt({n}) - {x}" / "{x}**2 - self.n")
    newton_raphson.formula = bisektion.formula  # Reuse the same formula for Newton-Raphson
    regula_falsi.formula = bisektion.formula  # Reuse the same formula for Regula Falsi

    newton_raphson.formula_derivative = "2 * {x}"  # Derivative of the formula (e.g. "2 * {x}")

    bisektion.accuracy = 0.001  # Desired accuracy (e.g. 1e-50 / 0.001)
    newton_raphson.accuracy = bisektion.accuracy  # Reuse the same accuracy for Newton-Raphson
    regula_falsi.accuracy = bisektion.accuracy  # Reuse the same accuracy for Regula Falsi

    bisektion.start_method()  # Start the method selection process

# write plotter class
# -maybe- write gui class for options