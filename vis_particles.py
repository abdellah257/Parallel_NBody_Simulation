
import matplotlib.pyplot as plt

''' A Python code to visualize the points at the end of a simulation'''

f = open("pdist.data", 'r')

X = []
Y = []
for line in f:
    line = line.split(",")
    x = float(line[0])
    y = float(line[1])
    X.append(int(x))
    Y.append(int(y))

print(X)
plt.scatter(X, Y, c='r', label='data')
plt.show()