import matplotlib.pyplot as plt

''' A Python code to draw the performence graphs'''

f = open("speed_seq.txt", 'r')
g = open("speed_par.txt", 'r')
X, Y = [], []
T, S = [], []
for line in f:
    line = line.split(",")
    x = float(line[0])
    y = float(line[2])
    X.append(int(x))
    Y.append(y)

for line in g:
    line = line.split(",")
    x = float(line[0])
    y = float(line[2])
    T.append(int(x))
    S.append(y)

plt.ylabel("Time in seconds")
plt.xlabel("N number of particles")
plt.plot(X, Y, c='b', label='Sequential')
plt.plot(T, S, c='r', label='Parallel')
plt.legend()
plt.show()