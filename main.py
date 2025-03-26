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
        self.controll_value:float = 0.0
        self.result:dict = {}
        self.formula:str = ""
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
        self.controll_value = self.fa * self.fb
        if self.controll_value < 0:
            self.controll = True
        else:
            self.controll = False
    
    def switch_interval(self):
        if self.fa * self.fc < 0:
            self.b = self.c
        else:
            self.a = self.c

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
        self.check_controll()
        
        print("="*28, "Iteration {iteration}".format(iteration = iteration), "="*28)
        print(f'a: {self.a}, b: {self.b}, c: {self.c}')
        print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
        print(f'controll: {self.controll}, controll value: {self.controll_value}')
        print("="*69)

        self.switch_interval()

        while abs(self.fc) > self.accuracy and self.controll:
            iteration += 1
            self.calculate_fa()
            self.calculate_fb()
            self.calculate_c()
            self.calculate_fc()
            self.check_controll()

            print("="*28, "Iteration {iteration}".format(iteration = iteration), "="*28)
            print(f'a: {self.a}, b: {self.b}, c: {self.c}')
            print(f'fa: {self.fa}, fb: {self.fb}, fc: {self.fc}')
            print(f'controll: {self.controll}, controll value: {self.controll_value}')
            print("="*70)

            self.switch_interval()

            if self.fc == 0:
                self.result["c:"] = self.c
                break
            elif self.fa == 0:
                self.result["a:"] = self.a
                break
            elif self.fb == 0:
                self.result["b:"] = self.b
                break
        
        if not self.controll:
            print("No solution found!")
        else:
            if self.fa * self.fc < 0:
                self.result["b:"] = self.b
                self.result["c:"] = self.c
            else:
                self.result["a:"] = self.a
                self.result["c:"] = self.c

        # fix result print
        print("="*31, "Result", "="*31)
        for key, value in self.result.items():
            print(f'{key} {value}')
        print("="*70)

if __name__ == "__main__":
    bisektion = Bisektion()
    bisektion.formula = "np.sqrt(self.n) - {x}"
    bisektion.accuracy = 0.001
    bisektion.main_loop()
