# 随机游走
import matplotlib.pyplot as plt
import random

random.seed(1)
position = 0
walk = [position]
steps = 200
for i in range(steps):
    step = 1 if random.randint(0, 1) else -1
    position += step
    walk.append(position)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(walk)
plt.title('One dimension Random Walk',size=20)
plt.xlabel('time')
plt.ylabel('location')
plt.savefig('1dim RW.png')
plt.show()
