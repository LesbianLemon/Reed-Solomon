class GaloisField: #Will represent finite fields (Z/p)[x]/f(x) (polynomials with coeficients from finite field Z/p (integers modulo p) modulo a prime polynomial f(x) of degree n) with size p^n written as GF(p^n)
  def __init__(self, p=2, n=8, alpha=2, prim_poly=285): #for Reed-Solomon codes we use GF(2^8) with polynomials [b7*x^7 + b6*x^6 + b5*x^5 + b4*x^4 + b3*x^3 + b2*x^2 + b1*x^1 + b0] represented as binary numbers [b7 b6 b5 b4 b3 b2 b1 b0]
    self.prime = p
    self.power = n
    self.alpha = alpha #primitive element from which all other elements can be derived (example: GF(8) = {0, 1, α, α^2, α^3, α^4, α^5, α^6})
    self.prim_poly = prim_poly #an irreducible primitive polynomial with which we build the Field

    self.size = p**n #amount of elements, self.size-1 is the cap value of our field
    self.expLUT = [-1]*(2*(self.size-1)) #lookup table for exponentiation of an element: expLUT[i] = α^i (we initialize a list of length 2*(self.size-1) to not have "IndexError: list index out of range" troubles in the future)
    self.logLUT = [-1]*self.size #lookup table for logarithm of an element: logLUT[a^i] = i (defined only for values above 0)

    a = 1
    for i in range(self.size-1): #we loop self.size-1 times, since we skip log(0) and since there are only self.size-1 different values in expLUT (there are self.size different elements in the field, but α^x; x∈GF(p^n) will never be 0, therefore they start repeating after self.size-1 iterations)
      if self.logLUT[a] != -1: #if this index has already been visited, we are not guaranteed unique values from α in this field, which is a necessity
        raise ValueError("the alpha and prim_poly arguments are not compatible")
      self.expLUT[i] = self.expLUT[i+(self.size-1)] = a #we build two identical tables in the same list with offset self.size-1 to avoid having to use modulo in the future
      self.logLUT[a] = i
      a = self.standard_mul(a, self.alpha) #αi = α(i-1)*α with α0 = 1 (1 -> α -> α^2 -> α^3...)

  def add(self, x, y):
    return x ^ y #assuming we are working with GF(2^n) (addition of polynomials works by adding coefficients with same degree (meaning no carry) and working in Z/2 is the same as modulo 2)

  def sub(self, x, y):
    return x ^ y #assuming we are working with GF(2^n) subtraction is the same as addition

  def standard_mul(self, x, y):
    result = 0
    while y > 0:
      if y & 1: #if y is odd
        result ^= x
      x <<= 1 #same as 2x
      y >>= 1 #same as y//2 (right shift drops the last bit)
      if x & self.size: #same as x > self.size-1 (self.size-1 indicates the maximum value that can exist in the field)
        x ^= self.prim_poly #if x gets larger than 255 (assuming GF(2^8)) it means the polynomial it represents is degree 8 or higher and must be reduced by our primitive polynomial
    return result

  def mul(self, x, y):
    if x == 0 or y == 0: #multiplication by 0
      return 0
    return self.expLUT[self.logLUT[x] + self.logLUT[y]] #x*y can be written as α^n*α^m = α^(n+m), where n and m are their log values

  def div(self, x, y):
    if y == 0: #division by 0
      raise ZeroDivisionError("cannot divide by zero")
    if x == 0: #zero divided by something
      return 0
    return self.expLUT[self.logLUT[x] - self.logLUT[y]] #x/y can be written as α^n/α^m = α^(n-m), where n and m are their log values