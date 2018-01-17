## Project 4 - Power Spectra and Signal Processing
## Physics 188B, Josh Samani
## Kevin Belleville, 14 June 2017

import matplotlib.pyplot as plt 
import numpy as np 

##################################################################################
## reading data and saving
##################################################################################

## third col of rays_monthly.txt
## corrected cosmic ray count rates

rays_data = []

with open("rays_monthly.txt", "r") as file:
	for line in file.readlines():
		stuff = line.split()
		rays_data.append(float(stuff[2]))

## second col of sunspots_monthly.txt
## monthly mean total sunspot number

sunspots_data = []

with open("sunspots_monthly.txt", "r") as file:
	for line in file.readlines():
		stuff = line.split()
		sunspots_data.append(float(stuff[1]))

## generate time array

N = 600
n = 1200
T = 600
time = np.linspace(0.0, 600, 600)

## generate freq array

k = np.arange(N)
freq = k/T
freq = freq[range(N//2)]


##################################################################################
## first assignment -- compute the power spectra of the data
##################################################################################

sunspots_fft = np.fft.rfft(sunspots_data)
sunspots_fft = sunspots_fft[range(N//2)]

power_sun = (abs(sunspots_fft)**2)*T*T/(N*N)

rays_fft = np.fft.rfft(rays_data)
rays_fft = rays_fft[range(N//2)]

power_rays = (abs(rays_fft)**2)*T*T/(N*N)


# plt.plot(time, sunspots_data, "r-")
# plt.show()
# plt.plot(freq, abs(sunspots_fft), "b-")
# plt.show()

# plt.plot(freq, power_sun, "y-")
# plt.show()

# plt.plot(time, rays_data, "r-")
# plt.show()
# plt.plot(freq, abs(rays_fft), "b-")
# plt.show()

# plt.plot(freq, power_rays, "y-")
# plt.show()


##################################################################################
## second assignment -- compute autocorrelation functions
##################################################################################

# normalize the data

def mean(data):
	""" Input: data set, list or array
		Returns: the mean of the data set, float
	"""
	sum = 0
	for d in data:
		sum += d
	return sum/(len(data))

def std(data, mu):
	""" Input: 	data set, list ;
				mu, float
		Returns: the standard deviation of the data set, float
	"""
	sum = 0
	for d in data:
		sum += (d - mu)**2
	sum = np.sqrt(sum)
	sum = sum/len(data)
	return sum

def normalize(d, mu, theta):
	""" Input: 	d; data point, float
				mu; mean, float
				theta; standard deviation, float
		Returns: normalized value, float
	"""
	new_value = d - mu
	new_value /= theta
	return new_value

# normalize sunspots
sunspots_normalized = []
sunspots_mean = mean(sunspots_data)
sunspots_std = std(sunspots_data, sunspots_mean)
for d in sunspots_data:
	_ = normalize(d, sunspots_mean, sunspots_std)
	sunspots_normalized.append(_)

# normalize rays
rays_normalized = []
rays_mean = mean(rays_data)
rays_std = std(rays_data, rays_mean)
for d in rays_data:
	_ = normalize(d, rays_mean, rays_std)
	rays_normalized.append(_)


# zero padding # stupid way of doing it but i already wrote it and it works

zero_padding = []
for i in range(N):
	zero_padding.append(0)

for zero in zero_padding:
	sunspots_normalized.append(zero)
	rays_normalized.append(zero)


##################################################################################
# correlation function implementation
##################################################################################

def l(n):
	if n < 600:
		return n
	else:
		return n-1200

def corr(x, y, n):
	""" Input:	x: first data set
				y: second data set
				n: point in the data set
		Returns: correlation data set: list"""
	first = np.fft.rfft(x)
	second = np.fft.rfft(y)
	# mult = [a * b for a,b in zip(first, second)]
	mult = first * second
	inv = np.fft.irfft(mult)
	ans = 1/(600-abs(l(n))) * inv
	return ans

# actually computing the auto correlation, sunspots
sunspots_autocorr = []
for i in range(1200):
	if i != 600:
		sunspots_autocorr.append(corr(sunspots_normalized, sunspots_normalized, i))

sunspots_autocorr = sunspots_autocorr[0]

# sunspots_autocorr = sunspots_autocorr[len(sunspots_autocorr)//2:]

# experimenting
# test = [0, 200, 400, 599, 601, 1000, 1188]
colors = ["r-", "b-", "g-", "y-", "m-", "c-", "k-"]
# for i, t in enumerate(test):
# 	_ = sunspots_autocorr[t]

x_values = np.linspace(0, 600, len(sunspots_autocorr))

# plt.plot(x_values, sunspots_autocorr, "c-")

# plt.show()

## numpy implementation of corr
# test = np.correlate(sunspots_normalized, sunspots_normalized, "full")
# test = test[test.size//2:]
# test_x = [i for i in range(len(test))]
# plt.plot(test_x, test)
# plt.show()

# actually computeing the auto correlation, rays
rays_autocorr = []
for i in range(1200):
	if i != 600:
		rays_autocorr.append(corr(rays_normalized, rays_normalized, i))

rays_autocorr = rays_autocorr[0]

rays_autocorr = rays_autocorr[len(rays_autocorr)//2:]

x_values_2 = np.linspace(0, 600, len(rays_autocorr))

# plt.plot(x_values_2, rays_autocorr, "m-")

# plt.show()

##################################################################################
# third assignment, using the cross-correlation function on sunspots AND rays
##################################################################################

cross_corr = []
for i in range(1200):
	if i != 600:
		cross_corr.append(corr(rays_normalized, sunspots_normalized, i))

cross_corr = cross_corr[0]
# cross_corr = cross_corr[len(cross_corr)//2:]

x_values_3 = np.linspace(-600, 600, len(cross_corr))

# plt.plot(x_values_3, cross_corr, "b-")
# plt.show()

##################################################################################
# fourth assignment, low-pass filter
# various cutoff filters
##################################################################################

sunspots_fft = np.fft.rfft(sunspots_data)
rays_fft = np.fft.rfft(rays_data)

def low_pass_filter(input_fft, cutoff):
	""" Parameters:	input_fft: the input fft, array or list
					cutoff: the cutoff frequency f_0, float
		Returns: a new array with a low pass filter applied to the input array """
	result = []
	if cutoff == 0:
		return np.fft.irfft(input_fft)
	for value in input_fft:
		if abs(value) <= cutoff:
			result.append(value)
		else:
			result.append(0)
	result = np.fft.irfft(result)
	return result

colors.append("-")
cutoff_freqs = [1000, 5000, 10000, 25000, 50000, 75000, 100000, 0]

t = np.linspace(0, 600, len(sunspots_data))


for i, f in enumerate(cutoff_freqs):
	plt.subplot(8,1,i+1)
	if i == 0:
		plt.title("Cutoff Frequency Comparisons")
	# plt.axis("off")
	if i == 4:
		plt.ylabel("Mean Sunspot Count")
	plt.plot(t, low_pass_filter(rays_fft, f), colors[i], label=str(f))
	# plt.ylim(-107, 299)
	plt.legend()
	if i != 7:
		plt.xticks([])


plt.xlabel("Time (months)")
plt.show()


# plt.plot(t, sunspots_data, "r-")
# plt.show()

# plt.plot(t, low_pass_filter(sunspots_fft, 250), "m-")
# plt.show()