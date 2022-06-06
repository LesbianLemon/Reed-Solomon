from galoisfield import GaloisField

class Polynomials:
  def __init__(self, field):
    self.field = field #the Galois Field from which the coefficients are taken from
    self.cap = field.size - 1 #cap value od coefficients

  def add(self, p, q):
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    for i in q:
      if i > self.cap:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    res = [0] * max(len(p), len(q)) #degree of new polynomial <= max(degree p, degree q)
    for i in range(len(p)):
      res[len(res) - len(p) + i] = p[i] #set coefficients of new polynomial as coefficients of polynomial p
    for i in range(len(q)):
      res[len(res) - len(q) + i] = self.field.add(res[len(res) - len(q) + i], q[i]) #add coefficients of new polynomial and coefficients of polynomial q
    return res

  def mul(self, p, q):
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    for i in q:
      if i > self.cap:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    res = [0] * (len(p) + len(q) - 1) #degree of new polynomial is degree p plus degree q
    for i in range(len(p)): #multiply every coefficient of p with every coefficient of q
      for j in range(len(q)):
        prod = self.field.mul(p[i], q[j]) #product of coefficient in Galois Field
        res[i + j] = self.field.add(res[i + j], prod) #add the product to the currect coefficient values based on degree
    return res

  def scalar(self, p, x):
    for i in p: #input must be constrained by Galois Field, most likely 0-255
      if i > self.cap:
        raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")
    if x > self.cap:
      raise ValueError("coefficients of given polynomials do not fit the charachteristics of the field")

    return [self.field.mul(x, coeff) for coeff in p]

# clss = Polynomials(GaloisField())
# print(clss.add([1, 255, 1], [2, 6]))
# print(clss.mul([1, 0, 2, 0], [1, 255, 3]))
# print(clss.scalar([1, 2, 3, 0, 1, 2, 3], 500))