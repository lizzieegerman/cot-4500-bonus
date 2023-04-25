from __init__ import *
import numpy as np

# Question 1:
np.set_printoptions(precision=5, suppress=True, linewidth=100)
A = np.array([[3, 1, 1], [1, 4, 1], [2, 3, 7]])
b = np.array([1, 3, 0])
toler = 1e-6
iter = 50

def norm(toler, x0, x):
    return (max(abs(x - x0))) / (max(abs(x0)) + toler)

def gauss_seidel(A, b, iter, toler):
    length = len(b)
    x = np.zeros(length)
    k = 1
    while (k <= iter):
        x0 = x.copy()
        for i in range(length):
            sum_1 = sum_2 = 0
            for j in range(i):
                sum_1 += (A[i][j] * x[j])
            
            for j in range(i + 1, length):
                sum_2 += (A[i][j] * (x0[j]))

            x[i] = (1 / A[i][i]) * (-sum_1 - sum_2 + b[i])

            if (norm(toler, x0, x) < toler):
                return k
        
        k += 1

    return k

print(gauss_seidel(A, b, iter, toler))
print()


# Question 2:
def jacobi(A, b, iter, toler):
    length = len(b)
    x = np.zeros(length)
    k = 1

    while (k <= iter):
        x0 = x.copy()
        for i in range(length):
            total = 0
            for j in range(length):
                if j != i:
                    total += (A[i][j] * x[j])

            x0[i] = (b[i] - total) / A[i][i]

        if (norm(toler, x0, x) < toler):
            return k
        x = x0

        k += 1
    return k

print(jacobi(A, b, iter, toler))
print()

# Question 3:
def f(x):
    return (x ** 3) - (x ** 2) + 2

def f_prime(x):
    return 3 * (x ** 2) - 2 * x

x = 0.5
toler = 1e-6
iteration = 100

def newton_raphson(x, toler, iter):
    k = 0
    while abs(f(x)) >= toler and k < iteration:
        quotient = f(x) / f_prime(x)
        x -= quotient
        k += 1
    print(k+1)  
    return k
newton_raphson(x, toler , iteration)
print()


# Question 4:
def applied_div_diff(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i + 2):
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue

            left: float= matrix[i][j-1]
            diagonal: float= matrix[i-1][j-1]
            numerator: float= left-diagonal
            denominator = matrix[i][0]-matrix[i-j+1][0]
            operation= numerator/denominator
            matrix[i][j]=operation
    return matrix      
        
def hermite_interpol():
    x_points = [0, 1, 2]
    y_points = [1, 2, 4]
    slopes = [1.06, 1.23, 1.55]
    num_of_points= len(x_points)
    matrix = np.zeros((2 * num_of_points, 2 * num_of_points))
    i = 0

    for x in range(0, num_of_points * 2, 2):
        matrix[x][0] = x_points[i]
        matrix[x+1][0] = x_points[i]
        i += 1


    i = 0    
    for x in range(0, num_of_points * 2, 2):
        matrix[x][1] = y_points[i]
        matrix[x+1][1] = y_points[i]
        i += 1


    i = 0   
    for x in range(1, num_of_points * 2, 2):
        matrix[x][2] = slopes[i]
        i += 1
    filled_matrix = applied_div_diff(matrix)
    print(filled_matrix)
hermite_interpol() 
print()


# Question 5:
def function(t, y):
    return y - t ** 3

def eulers_method():
    a = 0
    b = 3
    n = 100
    w = 0.5
    step = (b - a) / n
    for i in range(n):
        w = w + step / 2 * (function(a, w) + function(a + step, w + (step * function(a, w))))
        a += step
    return round(w, 5)

print("%.5f" % eulers_method())
print()