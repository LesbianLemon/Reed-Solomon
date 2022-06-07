from galoisfield import GaloisField
from polynomials import Polynomials

class ReedSolomonCodec:
  def __init__(self, enc_len):
    self.enc_len = enc_len

    self.field = GaloisField(2, 8, 2, 285) #field in which the encoding will take place
    self.polynomials = Polynomials(self.field) #holder class for operations with polynomials inside a Galois Field

  def generator_poly(self):
    gen = [1] #initial polynomial is just a constant (1) in case enc_len is equal to 0. gen_0 = 1
    for i in range(self.enc_len):
      gen = self.polynomials.mul(gen, [1, self.field.sub(0, self.field.pow(self.field.alpha, i))]) #gen_i = gen_(i-1) * (x - α^i) (for example: gen_4 = (x - 1)(x - α)(x - α^2)(x - α^3))
    return gen

  def encode(self):
    pass

# clss = ReedSolomonCodec(4)
# print(clss.generator_poly())