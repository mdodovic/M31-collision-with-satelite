import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from scipy import odr


def func(p, x):
    """
    Function that fit data
    Input parameters:
    p - array of parameters that needs to be founded
    x - variable
    Output parameter:
    analitycal funtion
    """
    a, b = p
    return a * x + b

# Model object with function
quad_model = odr.Model(func) 

# data and error
#real data and its error; arrays
x0 = [2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
y0 = [-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327, -0.886]
noise_x = [0.33] * 8
noise_y = [0.22] * 8
x0 = np.asarray(x0)
y0 = np.asarray(y0)
noise_x = np.asarray(noise_x)
noise_y = np.asarray(noise_y)
#Data with noise that fitted; array
y = y0 + noise_y
x = x0 + noise_x

x = np.asarray(x)
y = np.asarray(y)

# Create a RealData object - with data and errors
data = odr.RealData(x, y, sx=noise_x, sy=noise_y)

# Set up ODR with the model and data; beta0 are initial value of parameters - same number as in function
odr = odr.ODR(data, quad_model, beta0=[0., 1.])

# Run the regression.
out = odr.run()

# print fit parameters and 1-sigma estimates
popt = out.beta  # parameters
perr = out.sd_beta  # error of parameters
print("fit parameter 1-sigma error:")
for i in range(len(popt)):
    print(str(popt[i]) + ' +- ' + str(perr[i]))

# prepare confidence level curves
nstd = 5. # to draw 5-sigma intervals
popt_up = popt + nstd * perr
popt_dw = popt - nstd * perr

# parameters: down limit; upper limit; number of dots   
x_fit = np.linspace(min(x), max(x), 8) #parametters for function  
fit = func(popt, x_fit) #value of function for parameters 

# nstd-sigma field
fit_up = func(popt_up, x_fit)
fit_dw= func(popt_dw, x_fit)

#plot
fig, ax = plt.subplots(1)
rcParams['font.size']= 12
errorbar(x, y, yerr=noise_y, xerr=noise_x, hold=True, ecolor='k', fmt='none')#, label='data')
xlabel("x [\N{DEGREE SIGN}]")
ylabel("y [\N{DEGREE SIGN}]")
xlim(5,-0.5)
ylim(-5,0.5)
#title('fit with error on both axis', fontsize=18)
plot(x_fit, fit, 'r', lw=2, label='best fit curve')
#plot(x0, y0, 'k-', lw=2, label='True curve')
#ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='5-sigma interval')
legend(loc='upper left',fontsize=18)
show()
