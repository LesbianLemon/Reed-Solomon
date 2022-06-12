from typing import Iterable, Union
from galoisfield import GaloisField
from polynomials import Polynomials

class ReedSolomonCodec:
  def __init__(self, msg_len: int, enc_len: int, n: int=8, alpha: int=2, prim_poly: int=285) -> None:
    if msg_len <= enc_len:
      raise ValueError("message length must be bigger than length of parity bits")
    self.msg_len = msg_len
    self.enc_len = enc_len

    self.field = GaloisField(2, n, alpha, prim_poly) #field in which the encoding will take place
    self.polynomials = Polynomials(self.field) #holder class for operations with polynomials inside a Galois Field

  def generator_poly(self) -> list[int]:
    gen = [1] #initial polynomial is just a constant (1) in case enc_len is equal to 0. gen_0 = 1
    for i in range(self.enc_len):
      gen = self.polynomials.mul(gen, [1, self.field.sub(0, self.field.pow(self.field.alpha, i))]) #gen_i = gen_(i-1) * (x - α^i) (for example: gen_4 = (x - 1)(x - α)(x - α^2)(x - α^3))
    return gen

  def encode(self, msg: Iterable[int]) -> Iterable[int]:
    if len(msg) + self.enc_len > self.field.cap: #encoded message must be smaller than the cap
      raise ValueError(f"given message is too long to encode, cap is {self.field.cap}")

    generator = self.generator_poly() #generator polynomial with which we will divide to get a remainder
    padded = msg + [0]*(len(generator)-1) #pad the message to make room for remainder of polynomial division of msg with generator
    res, remainder = self.polynomials.monic_div(padded, generator) #only interested in the remainder as that is the encoding
    return msg + remainder #returning encoded message equal to res(x)*gen(x) in GF(2^n), where res(x) and gen(x) are polynomials for the result of division with generator and the generator polynomial itself

  def syndromes_poly(self, msg: Iterable[int]) -> list[int]:
    return [self.polynomials.eval(msg, self.field.pow(self.field.alpha, i)) for i in range(self.enc_len)][::-1] #evaluating the message with values used to create the generator polynomial, therefore if the recieved message is correct all coefficients should equal 0 (since a message without errors is just res(x)*gen(x))

  def sigma(self, pos: Iterable[int]) -> list[int]:
    sigma = [1] #initializing the sigma polynomial
    for i in pos:
      sigma = self.polynomials.mul(sigma, self.polynomials.sub([1], [self.field.pow(self.field.alpha, i), 0])) #formula for the sigma polynomial is as follows: sigma = sigma*(1 - alpha^i*x) for i in pos
    return sigma

  def omega(self, syndromes: Iterable[int], sigma: Iterable[int]) -> list[int]:
    return self.polynomials.monic_div(self.polynomials.mul(syndromes, sigma), [1] + [0]*(self.enc_len+1))[1] #omega(x) = (syndromes(x)*sigma(x)) mod x^(enc_len+1)

  def erasure_correction(self, msg: Iterable[int], pos: Iterable[int]) -> list[int]: #Forney algorithm - https://en.wikipedia.org/wiki/Forney_algorithm
    pos = [(len(msg)-1)-i for i in pos] #error positions must be reversed for the calculation of sigma(x)
    sigma = self.sigma(pos)
    omega = self.omega(self.syndromes_poly(msg), sigma)
    

clss = ReedSolomonCodec(70, 10)
print(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]))
# a = clss.polynomials.mul([clss.field.pow(2, 18), 1], [clss.field.pow(2, 15), 1])
# print(clss.polynomials.monic_div(clss.polynomials.mul(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]), a), [1] + [0]*(clss.enc_len+1)))