import numpy as np

class Bisektion():
    def __init__(self):
        self.a:float = 0.0
        self.b:float = 0.0
        self.c:float = 0.0
        self.n:float = 0.0
        self.fa:float = 0.0
        self.fb:float = 0.0
        self.fc:float = 0.0
        self.accuracy:float = 0.0
        self.formula:str = "np.sqrt(self.n) - {x}"
        self.controll:bool = False

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
                flag = False
            except ValueError:
                print("Please enter a valid number!")
                continue

    def calculate_fa(self):
        self.fa = eval(self.formula.format(x = self.a))

    def calculate_fb(self):
        self.fb = eval(self.formula.format(x = self.b))
    
    def calculate_c(self):
        self.c = (self.a + self.b) / 2
    
    def calculate_fc(self):
        self.fc = eval(self.formula.format(x = self.c))
    
    def check_controll(self):
        if self.fa * self.fc < 0:
            self.controll = True
        else:
            self.controll = False
    
    def switch_interval(self):
        ...

    def main_loop(self):
        ...

if __name__ == "__main__":
    bisektion = Bisektion()
    # move to main loop
    bisektion.get_n()
    bisektion.get_interval_a()
    bisektion.get_interval_b()
    bisektion.calculate_fa()