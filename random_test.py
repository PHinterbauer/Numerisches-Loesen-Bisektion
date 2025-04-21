class Base:
    def __init__(self, arg1: str = "default", arg2: int = 0, arg3: list = []):
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3

    def test(self):
        self.arg1 = "test"
        self.arg2 = 1
        self.arg3 = [1, 2, 3]

class Derived(Base):
    def __init__(self, base: Base, arg1: str = "der", arg2: int = 3, arg3: list = ["new"]):
        self.base = base
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3
    
    def __getattr__(self, attr):
        return getattr(self.base, attr)
    
class Derived2(Base):
    def __init__(self, base: Base):
        self.base = base
    
    def __getattr__(self, attr):
        return getattr(self.base, attr)

if __name__ == "__main__":
    base_instance = Base()
    derived_instance = Derived(base_instance)
    derived_instance2 = Derived2(derived_instance)
    print(base_instance.arg1, base_instance.arg2, base_instance.arg3)
    print(derived_instance.arg1, derived_instance.arg2, derived_instance.arg3)
    print(derived_instance2.arg1, derived_instance2.arg2, derived_instance2.arg3)
    # print("===========================================================================")
    # base_instance.test()
    # print(base_instance.arg1, base_instance.arg2, base_instance.arg3)
    # print(derived_instance.arg1, derived_instance.arg2, derived_instance.arg3)
    # print(derived_instance2.arg1, derived_instance2.arg2, derived_instance2.arg3)
    print("===========================================================================")
    derived_instance.test()
    print(base_instance.arg1, base_instance.arg2, base_instance.arg3)
    print(derived_instance.arg1, derived_instance.arg2, derived_instance.arg3)
    print(derived_instance2.arg1, derived_instance2.arg2, derived_instance2.arg3)

