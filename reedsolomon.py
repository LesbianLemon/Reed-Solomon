from typing import Iterable, Union
from galoisfield import GaloisField
from polynomials import Polynomials

class ReedSolomonCodec:
  def __init__(self, msg_len: int, enc_len: int, p: int=2, n: int=8, alpha: int=2, prim_poly: int=285) -> None:
    if msg_len <= enc_len:
      raise ValueError("message length must be bigger than length of parity bits")
    self.msg_len = msg_len
    self.enc_len = enc_len

    self.field = GaloisField(p, n, alpha, prim_poly) #field in which the encoding will take place
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
    syndromes = [0]*self.enc_len #initializing syndromes polynomial
    for i in range(0, self.enc_len):
      syndromes[i] = self.polynomials.eval(msg, self.field.pow(self.field.alpha, i)) #evaluating the recieved message with values used to create the generator polynomial. Since the encoded message (without any errors) in GF(2^n) equals res(x)*gen(x) all these evaluations should equal 0 as gen(x) will be 0 for these specific values
    return syndromes[::-1] #computer polynomial is in reverse order, highest term has to have highest power of alpha

clss = ReedSolomonCodec(70, 10)
# gen = clss.generator_poly()
# msg = [1, 2, 3, 4, 5]
# padded = [1, 2, 3, 4, 5, 0, 0, 0, 0]
# p, rem = clss.polynomials.monic_div(padded, gen)
# print(clss.polynomials.mul(p, gen), p, rem)
# print(clss.polynomials.add(clss.polynomials.mul(p, gen), rem))
# print(clss.encode(msg))
# print(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]))
# a = clss.polynomials.mul([clss.field.pow(2, 18), 1], [clss.field.pow(2, 15), 1])
# print(clss.polynomials.monic_div(clss.polynomials.mul(clss.syndromes_poly([0, 210, 117, 71, 118, 23, 50, 6, 39, 38, 150, 198, 198, 150, 112, 236, 188, 42, 144, 19, 107, 175, 239, 253, 75, 224]), a), [1] + [0]*(clss.enc_len+1)))