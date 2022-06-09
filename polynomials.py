from typing import Iterable

class Polynomials:
  def __init__(self, field):
    self.field = field #the Galois Field from which the coefficients are taken from
    self.cap = field.size - 1 #cap value od coefficients

  def add(self, p: Iterable[int], q: Iterable[int]) -> list[int]:
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    for i in q:
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    res = [0] * max(len(p), len(q)) #degree of new polynomial <= max(degree p, degree q)
    for i in range(len(p)):
      res[len(res) - len(p) + i] = p[i] #set coefficients of new polynomial as coefficients of polynomial p
    for i in range(len(q)):
      res[len(res) - len(q) + i] = self.field.add(res[len(res) - len(q) + i], q[i]) #add coefficients of new polynomial and coefficients of polynomial q
    return res

  def mul(self, p: Iterable[int], q: Iterable[int]) -> list[int]:
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    for i in q:
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    res = [0] * (len(p) + len(q) - 1) #degree of new polynomial is degree p plus degree q
    for i in range(len(p)): #multiply every coefficient of p with every coefficient of q
      for j in range(len(q)):
        prod = self.field.mul(p[i], q[j]) #product of coefficient in Galois Field
        res[i + j] = self.field.add(res[i + j], prod) #add the product to the currect coefficient values based on degree
    return res

  def scalar(self, p: Iterable[int], x: int) -> list[int]:
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    if x > self.cap or x < 0:
      raise ValueError("given scalar does not fit the charachteristics of the field")

    return [self.field.mul(x, coeff) for coeff in p]

  def monic_div(self, p: Iterable[int], q: Iterable[int]) -> list[int]: #expanded syntetic division with monic polynomials (expanded Horner's method) - https://en.wikipedia.org/wiki/Synthetic_division#Expanded_synthetic_division
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    for i in q:
      if i > self.cap or i < 0:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    if len(q) and q[0] != 1: #function only works with monic divisors (highest term coefficient must be 1)
      raise ValueError("given divisor is not monic")
    if len(p) < len(q):
      raise ValueError("given divisor is higher degree than dividend")

    res = p.copy()
    for i in range(len(p) - (len(q) - 1)): #expanded syntetic division has (degree p - (degree q - 1)) steps, since (degree q - 1) is the maximum degree of the remainder. The coefficients of higher degree parts can always be removed via division
      for j in range(1, len(q)): #skip first coefficient as the divisor is assumed to always be monic (first coefficient is always 1)
        res[i + j] = self.field.sub(res[i + j], self.field.mul(q[j], res[i])) #add negated j-th coefficient of q multiplied by i-th coeffcient of res to positions right of i
    separation = len(q) - 1
    return res[:-separation], res[-separation:]

# clss = Polynomials(GaloisField())
# print(clss.add([1, 255, 1], [2, 6]))
# print(clss.mul([1, 0, 2, 0], [1, 255, 3]))
# print(clss.scalar([1, 2, 3, 0, 1, 2, 3], 500))