from galoisfield import GaloisField

class ReedSolomonCodec:
  def __init__(self, enc_len):
    self.field = GaloisField(2, 8, 2, 285)

  def generator_poly(self, parity):
    pass

  def encode(self):
    pass