from typing import Iterable, Union
from galoisfield import GaloisField
from polynomials import Polynomials

class ReedSolomonCodec:
  """
  Create a Reed-Solomon code codec with a specified length of parity bits and Galois Field. When specifying a Galois Field the parameter 'n' creates a field GF(2^n), the parameters 'alpha' and 'prim_poly' written as a polynomials when in binary form (the number 11 or 1011 in binary refers to the polynomial [x^3 + x + 1]) change the characteristics of the field.

  Default Galois Field is GF(2^8) with a primitive element alpha as the polynomial [x] and a primitive polynomial prim_poly as [x^8 + x^4 + x^3 + x^2 + 1].

  Note: When changing 'alpha' and 'prim_poly' parameters make sure the polynomial represented by 'prim_poly' is irreducible and degree 'n', also confirm that the 'alpha' polynomial can truly generate all the elements in the field.
  """
  def __init__(self, msg_len: int, enc_len: int, n: int=8, alpha: int=2, prim_poly: int=285) -> None:
    if msg_len <= enc_len:
      raise ValueError("message length must be bigger than length of parity bits")
    self.msg_len = msg_len
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
      sigma = self.polynomials.mul(sigma, self.polynomials.sub([1], [self.field.pow(self.field.alpha, i), 0])) #formula for the sigma polynomial is as follows: sigma = sigma*(1 - alpha^i*x) for i in pos
    return sigma

  def omega(self, syndromes: Iterable[int], sigma: Iterable[int]) -> list[int]:
    """
    Return omega(x), a polynomial in list form (coefficients arranged from highest term to lowest) known as the error evaluator polynomial. It is calculated based on the syndromes polynomial, sigma(x) and amount of parity bits (parity bit amount is specified at initialization).
    """
    return self.polynomials.monic_div(self.polynomials.mul(syndromes, sigma), [1] + [0]*(self.enc_len+1))[1] #omega(x) = (syndromes(x)*sigma(x)) mod x^(enc_len+1)

  def __single_encode(self, msg: Iterable[int]) -> Iterable[int]:
    if len(msg) + self.enc_len > self.field.cap: #encoded message must be smaller than the cap
      raise ValueError(f"given message is too long to encode, cap is {self.field.cap}")

    padded = msg + [0]*self.enc_len #pad the message to make room for remainder of polynomial division of msg with generator
    remainder = self.polynomials.monic_div(padded, self.generator_poly())[1] #only interested in the remainder as that is the encoding
    return msg + remainder #returning encoded message equal to res(x)*gen(x) in GF(2^n), where res(x) and gen(x) are polynomials for the result of division with generator and the generator polynomial itself

  def encode(self, msg: Union[Iterable[int], str], return_str: bool=False) -> Union[Iterable[int], str]:
    str_marker = isinstance(msg, str)
    msg = [ord(ch) for ch in msg] if str_marker else msg
    encoded = self.__single_encode(msg)
    return "".join([chr(i) for i in encoded]) if return_str else encoded

  def erasure_correction(self, msg: Iterable[int], pos: Iterable[int]) -> list[int]: #Forney algorithm - https://en.wikipedia.org/wiki/Forney_algorithm
    """
    Return a polynomial containing the error factors at each term calculated with the usage of the recieved message and the positions of errors in the message via the Forney algorithm.
    
    When subtracted from the recieved message, the original message is recovered.
    """
    pos = [(len(msg)-1)-i for i in pos] #error positions must be reversed for the calculation of sigma(x)
    sigma = self.sigma(pos)
    if sigma.count(0) == len(sigma): #if message is without errors return polynomial with error factors 0
      return [0]*len(msg)
    omega = self.omega(self.syndromes_poly(msg), sigma)
    

# clss = ReedSolomonCodec(70, 4)
# print(clss.encode("hello world!"))
# print(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]))
# a = clss.polynomials.mul([clss.field.pow(2, 18), 1], [clss.field.pow(2, 15), 1])
# print(clss.polynomials.monic_div(clss.polynomials.mul(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]), a), [1] + [0]*(clss.enc_len+1)))