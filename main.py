import numpy as np
from plot import Plot

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
        self.result: dict = {}
        self.formula: str = ""
        self.controll: bool = False

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
                flag = False
            except ValueError:
                print("Please enter a valid number!")  # Handle invalid input
                continue

    def calculate_fa(self):
        """
        **calculate_fa**: Calculates the function value at a (fa).
        """
        self.fa = eval(self.formula.format(x=self.a))  # Evaluate formula at a

    def calculate_fb(self):
        """
        **calculate_fb**: Calculates the function value at b (fb).
        """
        self.fb = eval(self.formula.format(x=self.b))  # Evaluate formula at b

    def calculate_c(self):
        """
        **calculate_c**: Calculates the midpoint of the interval (c).
        """
        self.c = (self.a + self.b) / 2  # Midpoint formula

    def calculate_fc(self):
        """
        **calculate_fc**: Calculates the function value at c (fc).
        """
        self.fc = eval(self.formula.format(x=self.c))  # Evaluate formula at c

    def check_controll(self):
        """
        **check_controll**: Checks if the interval is valid by evaluating fa * fb.
        """
        self.controll_value = self.fa * self.fb  # Product of fa and fb
        self.controll = self.controll_value < 0  # Valid if product is negative

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
        while abs(self.fc) > self.accuracy and self.controll:
            iteration += 1
            self.calculate_fa()  # Update fa
            self.calculate_fb()  # Update fb
            self.calculate_c()  # Update c
            self.calculate_fc()  # Update fc
            self.check_controll()  # Recheck interval validity

            # Print iteration details
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
                self.result["c:"] = self.c  # Exact solution at c
                break
            elif self.fa == 0:
                self.result["a:"] = self.a  # Exact solution at a
                break
            elif self.fb == 0:
                self.result["b:"] = self.b  # Exact solution at b
                break

        # Handle cases where no exact solution is found
        if not self.result and self.controll:
            if self.fa * self.fc < 0:
                self.result["a:"] = self.a
                self.result["c:"] = self.c
            else:
                self.result["b:"] = self.b
                self.result["c:"] = self.c

        # Print the result
        if self.result:
            left_side_length, right_side_length, total_length = self.calc_separator(
                text_1=f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}',
                header=f' Result '
            )
            print("=" * left_side_length, "Result", "=" * right_side_length)
            for key, value in self.result.items():
                result_text = f'{key} {value}'
                left_padding, right_padding, _ = self.calc_separator(text_1="=" * total_length, header=result_text)
                print(" " * left_padding + result_text + " " * right_padding)
            print("=" * total_length)
        else:
            print("No solution found!")  # No valid solution

if __name__ == "__main__":
    bisektion = Bisektion()
    bisektion.formula = "np.sqrt(self.n) - {x}"  # Formula for root finding (e.g. "np.sqrt(self.n) - {x}" / "{x}**2 - self.n")
    bisektion.accuracy = 1e-50  # Desired accuracy (e.g. 1e-50 / 0.001)
    bisektion.main_loop()  # Start the bisection method
