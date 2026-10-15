prices = []

N = 60*60*24*365

for i in range(N):
    prices.append(100 + i/100)

print(prices[:5])

avg = 0.0

for p in prices:
    avg += p/len(prices)


import numpy as np

prices = 100 + np.arrange(N)/100
print(prices[:5])

