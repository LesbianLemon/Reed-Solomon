class GaloisField: #Will represent finite fields (Z/p)[x]/f(x) (polynomials with coeficients from finite field Z/p (integers modulo p) modulo a prime polynomial f(x) of degree n) with size p^n written as GF(p^n)
  """
  Create a Galois Field GF(p^n) with a primitive element 'alpha' and a primitive polynomial 'prim_poly'. The parameter 'p' refers to the prime with which we construct element polynomials of the field using coefficients from Z/p (integers modulo 'p'). The parameter 'n' is the power we raise the prime number to recieve the field size/cap (the highest number representing a polynomial). The parameters 'alpha' and 'prim_poly' are numbers refering to polynomials (like elements of GF) when expressed in number system 'p' (example: when p=2, the number 11 or 1011 in binary represents [x^3 + x + 1]). 'prim_poly' must be an irreducible polynomial of degree 'n' with we modulo to create the field. 'alpha' must be a primitive element of the field (meaning it can create all the other elements when raised to different powers).

  Default Galois Field is GF(2^8) with a primitive element alpha as the polynomial [x] and a primitive polynomial prim_poly as [x^8 + x^4 + x^3 + x^2 + 1].
  """
  def __init__(self, p: int=2, n: int=8, alpha: int=2, prim_poly: int=285) -> None: #for Reed-Solomon codes we use GF(2^8) with polynomials [b7*x^7 + b6*x^6 + b5*x^5 + b4*x^4 + b3*x^3 + b2*x^2 + b1*x^1 + b0] represented as binary numbers [b7 b6 b5 b4 b3 b2 b1 b0]
    if self.is_prime(p): #returns False for non-prime and values < 2
      if p > 2:
        raise NotImplementedError("current GaloisField class only supports GF(2^n)") #functionality still not added
      self.prime = p
    else:
      raise ValueError("enter a valid prime number")
    
    if n > 0:
      self.power = n
    else:
      raise ValueError("enter a valid power (positive integer)")

    self.alpha = alpha #primitive element from which all other elements can be derived (example: GF(8) = {0, 1, α, α^2, α^3, α^4, α^5, α^6})
    self.prim_poly = prim_poly #an irreducible primitive polynomial with which we build the Field

    self.size = p**n #amount of elements, self.size-1 is the cap value of our field
    self.cap = self.size - 1 #max element in field
    self.expLUT = [-1]*(2*(self.size-1)) #lookup table for exponentiation of an element: expLUT[i] = α^i (we initialize a list of length 2*(self.size-1) to not have "IndexError: list index out of range" troubles in the future)
    self.logLUT = [-1]*self.size #lookup table for logarithm of an element: logLUT[a^i] = i (defined only for values above 0)

    a = 1
    for i in range(self.size-1): #we loop self.size-1 times, since we skip log(0) and since there are only self.size-1 different values in expLUT (there are self.size different elements in the field, but α^x; x∈GF(p^n) will never be 0, therefore they start repeating after self.size-1 iterations)
      if self.logLUT[a] != -1: #if this index has already been visited, we are not guaranteed unique values from α in this field, which is a necessity
        raise ValueError("the alpha and prim_poly arguments are not compatible")
      self.expLUT[i] = self.expLUT[i+(self.size-1)] = a #we build two identical tables in the same list with offset self.size-1 to avoid having to use modulo in the future
      self.logLUT[a] = i
      a = self.standard_mul(a, self.alpha) #α_i = α_(i-1)*α with α_0 = 1 | (1 -> α -> α^2 -> α^3...)

  @staticmethod
  def is_prime(x: int) -> bool: #O(sqrt(n)) algorithm for prime checking
    """
    Return True if input x is a prime number, otherwise False.
    """
    if x < 2: #remove negative numbers, 0 and 1
      return False
    i = 2 #start at 2, since everything is divisible by 1
    while i*i <= x: #every non-prime number x has at least 2 divisors (besides 1 and itself), if one is found the number is not prime. Therefore we can simply check up until sqrt(x), since higher than that and the number cannot have 2 divisors
      if x % i == 0: #check if number is divisor
        return False #if any divisor is found, the number is not prime
      i += 1
    return True #no divisors found

  def add(self, x: int, y: int) -> int:
    """
    Return x+y in the Galois Field.

    Note: Currently only supported for GF(2^n).
    """
    if self.prime == 2:
      return x ^ y #assuming we are working with GF(2^n) (addition of polynomials works by adding coefficients with same degree (meaning no carry) and working in Z/2 is the same as modulo 2)
    raise NotImplementedError("the add method is currently only functional for GF(2^n)")

  def sub(self, x: int, y: int) -> int:
    """
    Return x-y in the Galois Field.

    Note: Currently only supported for GF(2^n)
    """
    if self.prime == 2:
      return x ^ y #assuming we are working with GF(2^n) subtraction is the same as addition
    raise NotImplementedError("the sub method is currently only functional for GF(2^n)")

  def standard_mul(self, x: int, y: int) -> int: #only for GF(2^n)
    """
    Return x*y in the Galois Field (the method is slower and only works for GF(2^n)).
    """
    if self.prime != 2:
      raise NotImplementedError(f"stardard_mul multiplication only works in GF(2^n) and not in GF({self.prime}^{self.power})")

    result = 0
    while y > 0:
      if y & 1: #if y is odd
        result = self.add(result, x)
      x <<= 1 #same as 2x
      y >>= 1 #same as y//2 (right shift drops the last bit)
      if x & self.size: #same as x > self.size-1 (self.size-1 indicates the maximum value that can exist in the field)
        x = self.add(x, self.prim_poly) #if x gets larger than 255 (assuming GF(2^8)) it means the polynomial it represents is degree 8 or higher and must be reduced by our primitive polynomial
    return result

  def mul(self, x: int, y: int) -> int:
    """
    Return x*y in the Galois Field.
    """
    if x == 0 or y == 0: #multiplication by 0
      return 0
    return self.expLUT[self.logLUT[x] + self.logLUT[y]] #x*y can be written as α^n*α^m = α^(n+m), where n and m are their log values

  def div(self, x: int, y: int) -> int:
    """
    Return x/y in the Galois Field.
    """
    if y == 0: #division by 0
      raise ZeroDivisionError("cannot divide by zero")
    if x == 0: #zero divided by something
      return 0
    return self.expLUT[self.logLUT[x] - self.logLUT[y]] #x/y can be written as α^n/α^m = α^(n-m), where n and m are their log value

  def pow(self, x: int, p: int) -> int:
    """
    Return x^p in the Galois Field.
    """
    return self.expLUT[(self.logLUT[x]*p) % self.cap]

  def inverse(self, x: int) -> int:
    """
    Return 1/x in the Galois Field.
    """
    return self.expLUT[-self.logLUT[x]] #x^(-1) can be written as α^(-n), with n being the log value

# clss = GaloisField()
# print(clss.mul(45, 216))
# print(clss.mul(54, 18))
# print(clss.inverse(2))