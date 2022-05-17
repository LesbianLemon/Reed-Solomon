class ByteMatrix:
  def __init__(self, rows, columns):
    self.rows = rows
    self.columns = columns

    self.data = []
    for i in range(rows):
      self.data.append(list([0]*columns))


  def set_pos(self, row, column, value):
    if row < 0 or row > self.rows:
      raise IndexError("row index out of range or negative")
    if column < 0 or column > self.columns:
      raise IndexError("column index out of range or negative")

    self.data[row][column] = value


  def set_data(self, data):
    rows = len(data)
    columns = len(data[0])

    if rows != self.rows:
      raise ValueError("row size does not match the one set in initialization")
    if columns != self.columns:
      raise ValueError("column size does not match the one set in initialization")

    for i in range(rows):
      if len(data[i]) != columns:
        raise ValueError("rows must be of equal length")
      for j in range(columns):
        self.set_pos(i, j, data[i][j])

  
  def get_pos(self, row, column):
    if row < 0 or row > self.rows:
      raise IndexError("row index out of range or negative")
    if column < 0 or column > self.columns:
      raise IndexError("column index out of range or negative")

    return self.data[row][column]

  
  def left_identity(self):
    result = ByteMatrix(self.rows, self.rows)
    for i in range(self.rows):
      result.set_pos(i, i, 1)
    return result


  def right_identity(self):
    result = ByteMatrix(self.columns, self.columns)
    for j in range(self.columns):
      result.set_pos(j, j, 1)
    return result


  def inverse(self):
    if self.rows != self.columns:
      raise ValueError("only square matricies have inverses")

    


  def __add__(self, right):
    if not isinstance(right, ByteMatrix):
      raise TypeError("unsupported addition with: " + str(type(right)))
    if self.rows != right.rows:
      raise ValueError("matricies must have same row size")
    if self.columns != right.columns:
      raise ValueError("matricies must have same column size")

    result = ByteMatrix(self.rows, self.columns)
    for i in range(self.rows):
      for j in range(self.columns):
        value = self.get_pos(i, j) + right.get_pos(i, j)
        result.set_pos(i, j, value)
    return result


  def __mul__(self, right):
    if not isinstance(right, ByteMatrix):
      raise TypeError("unsupported multiplication with: " + str(type(right)))
    if self.columns != right.rows:
      raise ValueError(f"right matrix must have row size (current: {right.rows}) equal to the column size (current: {self.columns}) of the left matrix")

    result = ByteMatrix(self.rows, right.columns)
    for i in range(self.rows):
      for j in range(right.columns):
        value = 0
        for k in range(self.columns):
          value += self.get_pos(i, k)*right.get_pos(k, j)
        result.set_pos(i, j, value)
    return result


  def __pow__(self, power):
    if self.columns != self.rows:
      raise ValueError(f"matrix must have row size (current: {self.rows}) equal to the column size (current: {self.columns})")
    
    temp = ByteMatrix(self.rows, self.columns)
    temp.set_data(self.data)
    result = ByteMatrix(self.rows, self.columns).left_identity()
    while power > 0:
      if power & 1:
        result *= temp
      temp *= temp
      power >>= 1
    return result
    

  def __str__(self):
    result = ""
    for i in range(self.rows):
      result += str(list(self.data[i])) + "\n"
    return result[:len(result)-1] #remove last newline


mat1 = ByteMatrix(4, 4)

mat1.set_data([
  [0, 1, 0, 1],
  [1, 0, 0, 1],
  [1, 0, 0, 0],
  [0, 1, 1, 0]
])

mat1 **= 2
print(mat1)

test = 0
for row in mat1.data:
  for el in row:
    test += el
print(test)