'''
Code adapted from the http://systems-sciences.uni-graz.at/etextbook/sw2/lyapunov.html , where you can find a very good description of the Lyapunov
exponent.

'''
import numpy as np
import matplotlib.pyplot as plt

lambdas = []

# define range of r 
r_values = np.linspace(0, 4, 100)

# Iterate through the different r values. 
for r in r_values:
    result = []
    x = 0.2

    # Iterate system 100 times and append the natural logarithms of the absolute value of F'(Xt)
    for i in range(100):
        result.append(np.log(np.abs((1 - r * x) * np.exp(r * (1 - x)))))
        x = x * np.exp(r * (1 - x))

    # Calculate lambda for each r value. Lambda is the arithmetic mean of the sumatory (sum(Xi)/2)
    lambdas.append(np.mean(result))  
    
fig = plt.figure(dpi=100, figsize=(6, 4))
ax1 = fig.add_subplot(1,1,1)

# Add lambda = 0 baseline
ax1.plot(r_values, [0]*len(r_values), 'g--')

# Plot Lyapunov exponent
ax1.plot(r_values, lambdas, 'b-', linewidth = 3, label = 'Lyapunov exponent (lambda)')

ax1.grid('on')
ax1.set_xlabel('r values')
ax1.legend(loc='best')
ax1.set_title('r values v.s. Lyapunov exponent')
plt.show