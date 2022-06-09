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

  def encode(self, msg: Union[Iterable[int], str]) -> Union[Iterable[int], str]:
    if len(msg) + self.enc_len > self.field.cap: #encoded message must be smaller than the cap
      raise ValueError(f"given message is too long to encode, cap is {self.field.cap}")

    if isinstance(msg, str): #if recieved input is in the form of a string, convert it to a polynomial with coefficients as ascii values
      tmp = msg
      msg = []
      for char in tmp:
        msg.append(ord(char)) #add ascii values as coefficients

    generator = self.generator_poly() #generator polynomial with which we will divide to get a remainder
    padded = msg + [0]*(len(generator)-1) #pad the message to make room for remainder of polynomial division of msg with generator
    res, remainder = self.polynomials.monic_div(padded, generator) #only interested in the remainder as that is the encoding
    
    encoded = msg + remainder
    # if self.field.size == 256:
    #   tmp = ""
    #   for i in encoded:
    #     tmp += chr(i)
    #   encoded = tmp
    return encoded

# clss = ReedSolomonCodec(70, 4)