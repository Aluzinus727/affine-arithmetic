from abc import ABC, abstractmethod

class Operable(ABC):
    @abstractmethod
    def to_intervals(self, intervals):
        pass

    @abstractmethod
    def get_range(self):
        pass

    @abstractmethod
    def __neg__(self):
        """-self"""
        pass

    @abstractmethod
    def __add__(self, other):
        """self + other"""
        pass
    
    @abstractmethod
    def __radd__(self, other):
        """other + self"""
        pass

    @abstractmethod
    def __sub__(self, other):
        """self - other"""
        pass
        
    @abstractmethod
    def __rsub__(self, other):
        """other - self"""
        pass

    @abstractmethod
    def __mul__(self, other):
        """self - other"""
        pass

    @abstractmethod
    def __rmul__(self, other):
        """self * other"""
        pass

    @abstractmethod
    def __rmul__(self, other):
        """other * self"""
        pass

    @abstractmethod
    def __truediv__(self, other):
        """self / other"""
        pass

    @abstractmethod
    def __rtruediv__(self, other):
        """other / self"""
        pass

    @abstractmethod
    def __pow__(self, other):
        pass
    
    @abstractmethod
    def cos(self, other):
        pass

    @abstractmethod
    def sin(self, other):
        pass

    @abstractmethod
    def sqrt(self, other):
        pass

    @abstractmethod
    def exp(self, other):
        pass

