from sage.libs.ntl import *
from sage.rings.polynomial.polynomial_integer_dense_ntl import *
import time

def set_ntl(element, modulus=None):
    # chooses between ZZX and ZZ_pX based on modulus given
    if modulus == 0 or modulus is None:
        return ntl.ZZX(element)
    else:
        return ntl.ZZ_pX(element, modulus)

class Poly:
    @classmethod
    def setup(cls, N=N, modulus=None):
        cls.N, cls.N2 = N, N*2
        # here we store the general modulus of the class, not to be confused with self.modulus
        cls.modulus = 0 if modulus is None else modulus
        cls.R = PolynomialRing(ZZ, 'x') # for printing
        cls.precomputations()
        
    @classmethod
    def precomputations(cls):
        cls.indices_auto5 = [ZZ((i * Zmod(cls.N2)(5)) % (cls.N2)) for i in range(cls.N)]
        cls.indices_auto5_poly = ntl.ZZ_pX(cls.indices_auto5, cls.N2)
              
    @classmethod
    def random(cls, modulus=None): # slow ~ 400ms, about 10x slower than NTL
        # maybe we use the NTL random function instead
        m = modulus if modulus is not None else cls.modulus
        ring = PolynomialRing(Zmod(m), 'x').quotient(x**cls.N + 1)
        return Poly(set_ntl(ring.random_element().list(), m), m)
    
    ## CREATION

    def __init__(self, coeffs, modulus=None):
        self.modulus = self.modulus if modulus is None else modulus
        
        if isinstance(coeffs, list):
            assert len(coeffs) <= self.N
            self.c = set_ntl(coeffs, modulus)
        elif isinstance(coeffs, Poly):
            self.c = coeffs.c
        else:    
            self.c = coeffs
            
    # ARITHMETIC OPERATORS
    
    def __add__(self, other):
        return Poly(self.c + other.c, self.modulus)
    
    def __radd__(self, other):
        return Poly(self.c + other.c, self.modulus)

    def __iadd__(self, other):
        self.c += other.c
        return self
    
    def __sub__(self, other):
        return Poly(self.c - other.c, self.modulus)
    
    def __rsub__(self, other):
        return Poly(other.c - self.c, self.modulus)
    
    def __isub__(self, other):
        self.c -= other.c
        return self

    def __neg__(self):
        return Poly(-self.c, self.modulus)
    
    ## MULTIPLICATION OPERATORS
    
    def __mul__(self, other): 
        if isinstance(other, Poly): # 90ms, if you square it's 60ms
            product = self.mod_quo(self.c * other.c)
        else: # integer multiplication, about 5ms
            product = self.c * set_ntl([other], self.modulus)
        return Poly(product, self.modulus)
    
    def __imul__(self, other):
        self.c *= other.c
        self.c = self.mod_quo(self.c)
        return self
    
    def __pow__(self, exponent):
        if exponent == 0:
            return Poly([1], self.modulus)
        elif exponent == 1:
            return self
        elif exponent % 2 == 0:
            result = self ** (exponent // 2)
            return result * result
        else:
            result = self ** ((exponent - 1) // 2)
            return self * result * result
        
    ## SCALING OPERATORS
    
    def rescale(self, other):
        return Poly(self.c._right_pshift(ntl.ZZ(other)), other)
    
    def __truediv__(self, other): # 5-6ms
        # in contrary to rescale, this does not scale down the modulus
        return self.rescale(other) % self.modulus        
    
    # MODULAR OPERATORS
    
    def __mod__(self, modulus): # fast, for the necessary cases 1-2ms
        Poly.modulus = modulus
        element_modulus = self.element_modulus()
        
        if modulus == 0: # this does convert to ZZX!!
            if element_modulus == 0:
                return self
            else: # slow takes 300ms
                return Poly(ntl.ZZX(self.c), 0)
            
        elif element_modulus == 0: # slow! 400ms
            return Poly(ntl.ZZ_pX(self.c, modulus), modulus)
        
        elif modulus == element_modulus:
            return self
        
        else:
            tmp = self.c.convert_to_modulus(ntl.ZZ_pContext(modulus))
            return Poly(tmp, modulus)
    
    def mod_quo(self, element, minus=True): # 1-2ms
        if minus: # mod X^N + 1
            return element.truncate(self.N) - element.left_shift(-self.N)
        else:
            return element.truncate(self.N) + element.left_shift(-self.N)
           
    ## SHIFTS AND AUTOMORPHISMS
        
    def __lshift__(self, n): # about 1-2ms
        temp = self.c.left_shift(n % self.N)
        return Poly(self.mod_quo(temp, minus=False), self.modulus)
    
    def __rshift__(self, n):
        return self << (self.N - n)
    
    def auto5(self): # 11 ms atm
        result = self.NTL_zero()
        for i in range(self.N):
            result[self.indices_auto5[i]] = self[i]
        return Poly(self.mod_quo(result), self.modulus)
    
    def auto(self, index): # 13 ms atm
        index = index % (self.N // 2)
        if index == 0: # index must be >= 1
            return self
        elif index == 1:
            return self.auto5()
        exponent = Zmod(self.N2)(5) ** (index - 1)
        # we use the NTL library to calculate the indices of the automorphism,
        # by storing the indices in a polynomial mod 2*N
        indices = self.indices_auto5_poly * ntl.ZZ_pX([exponent], self.N2)
        # same as auto5, but with a different permutation
        result = self.NTL_zero()
        for i in range(self.N):
            result[indices[i]] = self[i]
        return Poly(self.mod_quo(result), self.modulus)
    
    def auto_inverse(self): # 3ms atm
        length = self.c.degree()
        result = -self.c.reverse().left_shift(self.N - length)
        result[0] = -result[self.N]
        return Poly(result.truncate(self.N), self.modulus)        
    
    # ACCESSORS
    
    def __setitem__(self, key, value):
        self.c[key] = value
        
    def __getitem__(self, key):
        return self.c[key]
    
    def leading_coefficient(self):
        return self.c.leading_coefficient()
    
    ## CHECKS
    
    def is_zero(self):
        return self.c.is_zero()
    
    def is_one(self):
        return self.c.is_one()
    
    def is_monic(self):
        return self.c.is_monic()
    
    ## PRINTING AND REPRESENTATION
    
    def centered_list(self): # we perform the central reduction in [-Q//2, Q//2)
        if self.modulus == 0:
            return self.c.list()
        temp = set_ntl((self % self.modulus).c.list()).list()
        return [a if a < self.modulus//2 else a - self.modulus for a in temp]
    
    def norm(self): # slow, 400ms
        return max([abs(a) for a in self.R(self.centered_list())])
    
    def __repr__(self): # doesn't need to be fast
        return str(self.R(self.centered_list()))

    # OTHER METHODS
    
    def __copy__(self):
        return Poly(copy(self.c), self.modulus)
    
    def __eq__(self, other):
        assert isinstance(other, Poly), "Cannot compare with non-Poly object!"
        assert other.modulus == self.modulus, "Different moduli!"
        other.check_modulus()
        self.check_modulus()
        return self.c == other.c
    
    def NTL_zero(self): # gives a zero for the current ring
        tmp = copy(self.c)
        tmp.clear()
        return tmp
    
    def zero(self): 
        return Poly(self.NTL_zero(), self.element_modulus())
        
    def list(self):
        return self.c.list()
    
    def clear(self): # Resets this polynomial to zero, changes in place
        self.c.clear()
        
    def element_modulus(self):
        try: ## ZZ_pX case
            return self.c.get_modulus_context().modulus()
        except: ## ZZX case
            return 0
    
    def check_modulus(self): # Checks if the elements modulus is the same as the NTL modulus
        tmp = (self.modulus == self.element_modulus())
        if not tmp:
            print(f"Warning: Element modulus = {self.modulus} is different from the NTLs modulus {self.element_modulus()}!")
        return tmp
