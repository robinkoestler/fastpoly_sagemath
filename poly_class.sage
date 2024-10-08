from sage.libs.ntl import *
from sage.rings.polynomial.polynomial_integer_dense_ntl import *
from sage.libs.ntl.ntl_ZZ_p import ntl_ZZ_p_random_element
from sage.structure.element import is_Element
import time
import numpy as np

def set_ntl(element, modulus=None):
    if type(element) == np.ndarray: element = list(element)
    # chooses between ZZX and ZZ_pX based on modulus given
    if modulus == 0 or modulus is None:
        return ntl.ZZX(element)
    else: return ntl.ZZ_pX(element, modulus)

class Poly:
    @classmethod
    def setup(cls, N=N, modulus=None):
        cls.N, cls.N2 = N, N*2
        # here we store the general modulus of the class, not to be confused with self.modulus
        cls.modulus = 0 if modulus is None else modulus
        if cls.modulus & cls.modulus - 1 != 0:
            print(f"Warning: Modulus {cls.modulus} should be a power of 2!")
        cls.R = PolynomialRing(ZZ, 'x') # for printing
        cls.precomputations()
        
    @classmethod
    def precomputations(cls):
        cls.indices_auto5 = [ZZ((i * Zmod(cls.N2)(5)) % (cls.N2)) for i in range(cls.N)]
        cls.indices_auto5_poly = ntl.ZZ_pX(cls.indices_auto5, cls.N2)
        cls.sum_of_monomials = Poly(ntl.ZZ_pX([1] * cls.N, cls.modulus))

    @classmethod
    def random(cls, modulus): # ˜40ms
        temp = set_ntl([0], modulus=modulus)
        for i in range(cls.N):
            temp[i] = ntl_ZZ_p_random_element(modulus)
        return Poly(temp, modulus)

    @classmethod
    def random2(cls, modulus): # 400ms, deprecated, for stupid reasons this is slower
        assert modulus & modulus - 1 == 0, "Modulus must be a power of 2!"
        int64 = 2**63
        quo = ZZ(floor(log(modulus // 2, int64)))
        int64_quo = int64**quo
        rest = modulus // int64_quo // 2
        result = Poly(set_ntl([0], modulus=modulus))
        for i in range(quo):
            a = np.random.randint(-int64, int64, size=cls.N)
            a = set_ntl(list(a), modulus)
            result += Poly(a, modulus) * (int64**i)
        if rest > 1:
            a = np.random.randint(-rest, rest, size=cls.N)
            a = set_ntl(list(a), modulus)
            result += Poly(a, modulus) * int64_quo
        return result
            
    
    ## CREATION

    def __init__(self, coeffs, modulus=None):
        assert hasattr(self, 'N'), "You need to call setup() first!"
        self.modulus = self.modulus if modulus is None else modulus
        if self.modulus & self.modulus - 1 != 0:
            print(f"Warning: Modulus {self.modulus} should be a power of 2!")
        if is_Element(coeffs):
            print("(Poly): Warning: You are passing a SageMath object as coefficients!")
            coeffs = coeffs.list()
        if isinstance(coeffs, list):
            self.c = set_ntl(coeffs, self.modulus)
        elif isinstance(coeffs, np.ndarray):
            self.c = set_ntl(np.round(coeffs).astype(int), self.modulus)
        elif isinstance(coeffs, Poly):
            self.c = coeffs.c
        else:
            self.c = coeffs
        assert self.c.degree() < self.N, f"Degree of polynomial is too high: {self.c.degree()} > {self.N} (setup degree)!"
            
    # ARITHMETIC OPERATORS
    
    def __add__(self, other):
        return Poly(self.c + other.c, self.modulus)
    
    def __radd__(self, other): return self + other

    def __iadd__(self, other):
        self.c += other.c
        return self
    
    def __sub__(self, other):
        return Poly(self.c - other.c, self.modulus)
    
    def __rsub__(self, other): return -self + other
    
    def __isub__(self, other):
        self.c -= other.c
        return self

    def __neg__(self): return Poly(-self.c, self.modulus)
    
    ## MULTIPLICATION OPERATORS
    
    def __mul__(self, other): 
        if isinstance(other, Poly): # 90ms, if you square it's 60ms
            product = self.mod_quo(self.c * other.c)
        else: # integer multiplication, about 5ms
            product = self.c * set_ntl([other], self.modulus)
        return Poly(product, self.modulus)
    
    def __rmul__(self, other):
        if isinstance(other, Poly):
            return other * self
        else:
            product = self.c * set_ntl([other], self.modulus)
        return Poly(product, self.modulus)
    
    def __imul__(self, other):
        self.c *= other.c
        self.c = self.mod_quo(self.c)
        return self
    
    def __pow__(self, exponent):
        if exponent == 0: return Poly([1], self.modulus)
        elif exponent == 1: return self
        elif exponent % 2 == 0:
            result = self ** (exponent // 2)
            return result * result
        else:
            result = self ** ((exponent - 1) // 2)
            return self * result * result
        
    ## SCALING OPERATORS
    
    def center_div(self, other):
        assert other != 0, "Division by zero!"
        if other == 1:
            return self
        if self.modulus == 0:
            for i in range(self.N):
                self[i] = ZZ(self[i]) // other
        else: # we need to center appropriately first, see .center() method
            half_modulus = self.modulus // 2 - 1
            for i in range(self.N):
                value = ZZ(self.c[i])
                if value > half_modulus:
                    self.c[i] = (value - self.modulus) // other
                else:
                    self.c[i] = value // other
        return self
        
    
    def __truediv__(self, other): # 5-6ms, deprecate
        if other == 1:
            return self
        # in contrary to scale(), this does not scale down the modulus, and it does not reduce
        # but it does not round as expected.
        assert self.modulus % other == 0, "Modulus must be divisible by the scaling factor!"
        return Poly(self.c._right_pshift(ntl.ZZ(other)), self.modulus)
    
    def scale(self, other, newmod=False, center=False): # 7ms
        # this scales and rounds correctly, that is centralized!
        if other == 1:
            return self
        assert self.modulus % other == 0, "Modulus must be divisible by the scaling factor!"
        assert self.modulus & self.modulus - 1 == 0, "Modulus must be a power of 2!"
        assert other & other - 1 == 0, "Scaling factor must be a power of 2!"
        if center: # We run into a wrong algorithm below, if the coefficients are negative within the modulus.
            self = self.center()
        shift = (self.sum_of_monomials * (other // 2)) % self.modulus # for the rounding
        result = (self + shift).c._right_pshift(ntl.ZZ(other))
        #assert not (center and newmod), "Both should not be used together"
        if newmod:
            quotient = self.modulus // other
            if center:
                return Poly(result, quotient).center()
            return Poly(result, quotient) % quotient
        if center:
            return Poly(result, self.modulus)
        return Poly(result, self.modulus) % self.modulus
    
    def newscale(self, other, newmod=False):
        new = self % 0
        for i in range(self.N):
            new[i] = new[i]._integer_() // other
        return new

    def newscale2(self, other, newmod=True):
        new = self.__copy__()
        new = new.center()
        for i in range(self.N):
            new[i] = self[i]._integer_() // other
        if newmod:
            new_modulus = self.modulus // other
            return Poly(new % new_modulus, new_modulus)
        else:
            return Poly(new % self.modulus, self.modulus)

    # MODULAR OPERATORS
    
    def __mod__(self, modulus): # fast, for the necessary cases 2-3ms
        # Avoid calling % p if the modulus is 0!
        
        Poly.modulus = modulus
        element_modulus = self.element_modulus()
        
        if modulus == 0: # this does convert to ZZX!!
            if element_modulus == 0: return self
            else: # ~ 8ms
                zero = ntl.ZZX([0])
                zero.preallocate_space(self.N)
                for i in range(self.N):
                    zero[i] = self[i].lift_centered()
                return Poly(zero, 0)
            
        elif element_modulus == 0: # slow! 333ms
            return Poly(set_ntl(self.c.list(), modulus), modulus)
        
        elif modulus == element_modulus: return self
        
        else: # 2-3ms
            tmp = self.c.convert_to_modulus(ntl.ZZ_pContext(modulus))
            return Poly(tmp, modulus)
        
    def center(self): # For a given modulus, this changes the polynomial to its centered representation, ˜24ms
        assert self.modulus > 0, "Modulus 0 is not allowed"
        half_modulus = self.modulus // 2 - 1
        for i in range(self.N):
            if self.c[i]._integer_() > half_modulus:
                self.c[i] -= self.modulus
        return self
                
    
    def mod_quo(self, element, minus=True): # mod X^N+1, 1-2ms
        if minus: return element.truncate(self.N) - element.left_shift(-self.N)
        else: return element.truncate(self.N) + element.left_shift(-self.N)
           
    ## SHIFTS
        
    def __lshift__(self, n, monom=False): # about 1-2ms
        # this rotates the coefficients
        temp = self.c.left_shift(n % self.N)
        return Poly(self.mod_quo(temp, minus=monom), self.modulus)
    
    def __rshift__(self, n): return self << (self.N - n)
    
    def monomial_shift(self, n): # 1-2ms
        return self.__lshift__(n, monom=True)

    def left_shift(self, n): # we shift the coefficients only!
        return Poly(self.c.left_shift(n).truncate(self.N), self.modulus)

    def right_shift(self, n): # we shift the coefficients only!
        return Poly(self.c.right_shift(n), self.modulus)

    ## AUTOMORPHISMS
    
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

    def imaginary_part(self):
        return (self.auto(-1) - self).monomial_shift(self.N // 2)

    # TRACE AND NORM

    def trace(self, index = None):
        if self.N != 1 and self.N != 2:
            assert self.N > 2, "Trace is only defined for N > 2!"
            assert self.N & self.N - 1 == 0, "N must be a power of 2!"
            if index is None:
                index = log(self.N // 2, 2)
            auto_index = self.N // 4
            for i in range(index):
                self = self + self.auto(auto_index)
                auto_index //= 2
        return self
    
    def prod_norm(self, index = None): # algebraic norm
        if self.N != 1 and self.N != 2:
            assert self.N > 2, "Norm is only defined for N > 2!"
            assert self.N & self.N - 1 == 0, "N must be a power of 2!"
            if index is None:
                index = log(self.N // 2, 2)
            auto_index = self.N // 4
            for i in range(index):
                self = self * self.auto(auto_index)
                auto_index //= 2
        return self
    
    # ACCESSORS

    def __len__(self): return self.N
    def __setitem__(self, key, value): self.c[key] = value
    def __getitem__(self, key): return self.c[key]
    def leading_coefficient(self): return self.c.leading_coefficient()
    
    ## CHECKS
    
    def is_zero(self): return self.c.is_zero()
    def is_one(self): return self.c.is_one()
    def is_monic(self): return self.c.is_monic()
    
    ## PRINTING AND REPRESENTATION

    def list(self, full=False):
        if full:
            l = self.c.list()
            return l + [0] * (self.N - len(l))
        return self.c.list()
    
    def centered_list(self, full=False): # ~ 20ms
        # we perform the central reduction in (-Q//2, Q//2]
        if self.modulus == 0:
            return self.list(full=full)
        result = [0] * self.N
        for i in range(self.N):
            result[i] = self.c[i].lift_centered()
        return result
    
    def norm(self): # ~30ms
        return max([abs(ZZ(a)) for a in self.centered_list()])
    
    def __repr__(self): # doesn't need to be fast
        return str(self.R(self.centered_list()))

    # OTHER METHODS
    
    def __copy__(self): # doesn't really work in Sage??
        return Poly(self.c.__copy__(), self.modulus)
    
    def __eq__(self, other):
        assert isinstance(other, Poly), "Cannot compare with non-Poly object!"
        assert other.modulus == self.modulus, "Different moduli!"
        other.check_modulus()
        self.check_modulus()
        return self.c == other.c
    
    def NTL_zero(self): # gives a zero for the current ring
        tmp = self.c.__copy__()
        tmp.clear()
        return tmp
    
    def zero(self): 
        return Poly(self.NTL_zero(), self.element_modulus())
    
    def const(self, value):
        m = self.element_modulus()
        return Poly(set_ntl([value], m), m)
    
    def monomial(self, index):
        return self.const(1) << index
    
    def reverse(self):
        result = self.c.reverse().left_shift(self.N - self.c.degree() - 1)
        return Poly(result.truncate(self.N), self.modulus)

    def to_int(self):
        for i in range(self.N):
            self[i] = self[i]._integer_()
        return self
    
    def lift(self):
        if self.modulus == 0:
            return [ZZ(i) for i in self.list(full=True)]
        else:
            print(".lift(): Not implemented!")
    
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
