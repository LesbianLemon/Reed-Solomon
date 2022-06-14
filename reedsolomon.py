from typing import Iterable, Union
from galoisfield import GaloisField
from polynomials import Polynomials

class ReedSolomonCodec:
  """
  Create a Reed-Solomon code codec with a specified length of parity bits and Galois Field. When specifying a Galois Field the parameter 'n' creates a field GF(2^n), the parameters 'alpha' and 'prim_poly' written as a polynomials when in binary form (the number 11 or 1011 in binary refers to the polynomial [x^3 + x + 1]) change the characteristics of the field.

  Default Galois Field is GF(2^8) with a primitive element alpha as the polynomial [x] and a primitive polynomial prim_poly as [x^8 + x^4 + x^3 + x^2 + 1].

  Note: When changing 'alpha' and 'prim_poly' parameters make sure the polynomial represented by 'prim_poly' is irreducible and degree 'n', also confirm that the 'alpha' polynomial can truly generate all the elements in the field.
  """
  def __init__(self, enc_len: int, n: int=8, alpha: int=2, prim_poly: int=285) -> None:
    self.enc_len = enc_len

    self.field = GaloisField(2, n, alpha, prim_poly) #field in which the encoding will take place
    self.polynomials = Polynomials(self.field) #holder class for operations with polynomials inside a Galois Field

  def generator_poly(self) -> list[int]:
    """
    Return the generator polynomial in list form (coefficients arranged from highest term to lowest) based on the amount of parity bits (parity bit amount is specified at initialization).

    A codec with parity bit amount n will have a generator polynomial of degree n (length of list will be n+1, since we also have the constant at the end).
    """
    gen = [1] #initial polynomial is just a constant (1) in case enc_len is equal to 0. gen_0 = 1
    for i in range(self.enc_len):
      gen = self.polynomials.mul(gen, [1, self.field.sub(0, self.field.pow(self.field.alpha, i))]) #gen_i = gen_(i-1) * (x - α^i) (for example: gen_4 = (x - 1)(x - α)(x - α^2)(x - α^3))
    return gen

  def syndromes_poly(self, msg: Iterable[int]) -> list[int]:
    """
    Return the syndromes polynomial in list form (coefficients arranged from highest term to lowest) based on the inputted message and amount of parity bits (parity bit amount is specified at initialization).

    If the inputted message is without errors, the returned polynomial will have all coefficients equal to 0.
    """
    return [self.polynomials.eval(msg, self.field.pow(self.field.alpha, i)) for i in range(self.enc_len)][::-1] #evaluating the message with values used to create the generator polynomial, therefore if the recieved message is correct all coefficients should equal 0 (since a message without errors is just res(x)*gen(x))

  def sigma(self, pos: Iterable[int]) -> list[int]:
    """
    Return sigma(x), a polynomial in list form (coefficients arranged from highest term to lowest) known as the error locator polynomial. It is calculated based on the inputted list of error positions.

    Note: The positions inputted follow a reversed 0-index notation (the last element is at position 0 and the first is the highest position in the message).
    """
    sigma = [1] #initializing the sigma polynomial
    for i in pos:
      sigma = self.polynomials.mul(sigma, self.polynomials.add([1], [self.field.pow(self.field.alpha, i), 0])) #formula for the sigma polynomial is as follows: sigma = sigma*(1 - alpha^i*x) for i in pos
    return sigma

  def omega(self, syndromes: Iterable[int], sigma: Iterable[int]) -> list[int]:
    """
    Return omega(x), a polynomial in list form (coefficients arranged from highest term to lowest) known as the error evaluator polynomial. It is calculated based on the syndromes polynomial, sigma(x) and amount of parity bits (parity bit amount is specified at initialization).
    """
    return self.polynomials.monic_div(self.polynomials.mul(syndromes, sigma), [1] + [0]*self.enc_len)[1] #omega(x) = (syndromes(x)*sigma(x)) mod x^(enc_len+1)

  def __single_encode(self, msg: Iterable[int]) -> Iterable[int]:
    """
    Return a polynomial in list form (coefficients arranged from highest term to lowest) representing the encoded message without any errors.

    Note: This method can recieve the message input only in polynomial form and can encode only messages that will not surpass the field cap after encoding.
    """
    if len(msg) + self.enc_len > self.field.cap: #encoded message must be smaller than the cap
      raise ValueError(f"given message is too long to encode, cap is {self.field.cap}")

    padded = msg + [0]*self.enc_len #pad the message to make room for remainder of polynomial division of msg with generator
    remainder = self.polynomials.monic_div(padded, self.generator_poly())[1] #only interested in the remainder as that is the encoding
    return msg + remainder #returning encoded message equal to res(x)*gen(x) in GF(2^n), where res(x) and gen(x) are polynomials for the result of division with generator and the generator polynomial itself

  def encode(self, msg: Union[Iterable[int], str], return_str: bool=False) -> Union[Iterable[int], str]:
    """
    Return a polynomial in list form (coefficients arranged from highest term to lowest) representing the encoded message without any errors.

    By setting the 'return_str' flag to True, the returned message will not be in polynomial form, but rather a string with the coefficients replaced by ASCII characters.
    """
    encoded = []
    slicing_idx = list(range(0, len(msg), (self.field.cap - self.enc_len))) + [len(msg)] #indexes used for slicing the message to correct lengths

    for i, j in zip(slicing_idx[:-1], slicing_idx[1:]): #iterate through pairs of values of slicing_idx with one shift (f.e. [(0, 1), (1, 2), (2, 3), (3, None)])
      slice = [ord(ch) for ch in msg[i:j]] if isinstance(msg, str) else msg[i:j] #change from str to polynomial for each slice
      encoded += self.__single_encode(slice)
    return "".join([chr(i) for i in encoded]) if return_str else encoded

  def error_poly(self, msg: Iterable[int], pos: Iterable[int]) -> list[int]: #Forney algorithm - https://en.wikipedia.org/wiki/Forney_algorithm
    """
    Return a polynomial in list form (coefficients arranged from highest term to lowest) containing the error factors at each term calculated with the usage of the recieved message and the positions of errors in the message via the Forney algorithm.
    
    When subtracted from the recieved message, the original message is recovered.
    """
    reversed_pos = [(len(msg)-1)-i for i in pos] #error positions must be reversed so that they match the degree of the term
    sigma = self.sigma(reversed_pos)
    if sigma.count(0) == len(sigma): #if message is without errors return polynomial with error factors 0
      return [0]*len(msg)
    omega = self.omega(self.syndromes_poly(msg), sigma)

    roots = []
    for i in reversed_pos:
      roots.append(self.field.inverse(self.field.pow(self.field.alpha, self.field.cap - i))) #α^(-(cap - i)) is the same as α^i, but it better describes what is going on. We are looking for the inverse of α raised to the position in the message

    error_poly = [0]*len(msg)
    for i in range(len(roots)):
      X = roots[i]
      X_inverse = self.field.inverse(X)
      
      denominator = list(sigma) #copy of sigma(x)
      denominator = [coef if i%2 != 0 else 0 for i, coef in enumerate(sigma[::-1])][::-1]
      denominator = self.polynomials.eval(self.polynomials.monic_div(denominator, [1, 0])[0], X_inverse)

      numerator = self.field.mul(self.polynomials.eval(omega, X_inverse), X)

      error_poly[pos[i]] = self.field.div(numerator, denominator)
    return error_poly

  def decode_erasures(self, msg: Union[Iterable[int], str], pos: Iterable[int], return_str: bool=False) -> Union[list[int], str]:
    """
    Return a polynomial in list form (coefficients arranged from highest term to lowest) representing the corrected message with the provided erasure positions.

    By setting the 'return_str' flag to True, the returned message will not be in polynomial form, but rather a string with the coefficients replaced by ASCII characters.
    """
    if len(pos) > self.enc_len:
      raise ValueError(f"Reed-Solomon codes can only correct up to amount of parity bits errors, currently your cap is {self.enc_len}")

    decoded = []
    slicing_idx = list(range(0, len(msg), self.field.cap)) + [len(msg)] #indexes used for slicing the message to correct lengths

    for i, j in zip(slicing_idx[:-1], slicing_idx[1:]): #iterate through pairs of values of slicing_idx with one shift (f.e. [(0, 1), (1, 2), (2, 3), (3, 4)...])
      if len([err for err in pos if err >= i and err <= j]) > self.enc_len:
        raise ValueError(f"too many erasures between index {i} and index {j}. This Codec can only correct up to {self.enc_len} erasures every {self.cap} characters, to fix this increase the cap or provide less erasures")

      slice = [ord(ch) for ch in msg[i:j]] if isinstance(msg, str) else msg[i:j]
      decoded += self.polynomials.add(slice, self.error_poly(slice, pos))
    return "".join([chr(i) for i in decoded]) if return_str else decoded

  def erasure_sim(self, msg: Iterable[int], pos: Iterable[int]) -> list[int]:
    """
    Return a polynomial in list form (coefficients arranged from highest term to lowest) with coefficients at provided erasure positions set to 0.
    """
    if max(pos) >= len(msg):
      raise ValueError("erasure position indexes out of range of the message")
    
    res = list(msg)
    for i in pos:
      res[i] = 0
    return res