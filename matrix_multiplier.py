import numpy as np

def naive_multiply(matrix1, matrix2):
    rows1, cols1 = matrix1.shape
    rows2, cols2 = matrix2.shape
    if cols1 != rows2:
        print("Matricies must have the property: cols1=rows2")
        print("Aborting")
        return
    element = 0
    new_matrix = np.empty([rows1, cols2])
    for i in range(cols1):
        for j in range(rows2):
            for k in range(cols2):
                element += matrix1[i, k]*matrix2[k, j]
        print(element)
        new_matrix[i, j] = element
        element = 0
    return new_matrix

if __name__ == "__main__":
    matrix1 = np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
    matrix2 = matrix1
    print(matrix1)
    print(matrix1[0, 1])
    print(matrix1[1, 2])
    #print(naive_multiply(matrix1, matrix2))
    #print(np.matrix(matrix1)*np.matrix(matrix2))
