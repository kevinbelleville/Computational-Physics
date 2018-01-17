import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

##########################################################
## global paraemters for the model
##########################################################

nt = 30                        # num of temp points
N = 50                         # size of lattice
equil_steps = 1000            # num of equil steps
mc_steps = 1000              # num of mc steps

temp = 1.0


##########################################################
## the model code
##########################################################

def initial_state(N):
    ''' gen random spin config for lattice size N x N with values -1 or 1 '''
    state = 2*np.random.randint(2, size=(N, N))-1
    return state

def mc_moves(config, N, temp):
    ''' execute mc moves using metropolis alg '''
    for i in range(N):
        for j in range(N):
            ## random site i in spot [a, b] that we assign to spin s
            a = np.random.randint(0, N)
            b = np.random.randint(0, N)
            s = config[a, b]
            ## calculating delta E = 2 * s * (top, bot, left, right, H=0)
            ## including mod N if a or b are at the boundaries
            ## they lap over to the other side
            neighbors = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
            energy_change = 2*s*neighbors
            ## if energy change is <= 0, accept
            if energy_change <= 0:
                s *= -1
            ## elif accept with prob A = e^(-energy_change/T)
            elif np.random.rand() < np.exp(-energy_change/temp):
                s *= -1
            ## apply change to the spot
            config[a, b] = s
    return config

##########################################################
## calculating quantitites
##########################################################

## quantities can all be calculated with energy and mag

def calc_energy(config):
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4.


def calc_mag(config):
    ## simply sum all values
    return np.sum(config)


Temp = np.linspace(1, 4, nt)
Energy = np.zeros(nt)
Magnetization = np.zeros(nt)
SpecificHeat = np.zeros(nt)
Susceptibility = np.zeros(nt)


##########################################################
## running the code / plotting
##########################################################

## gen model
model2 = initial_state(N)
print(model2)

count = 0
for i in range(equil_steps):
    mc_moves(model2, N, temp)
    print(count)
    count += 1

## plotting stuff
cmap = mpl.colors.ListedColormap(['yellow', 'blue'])
bounds = [-N, -N, N, N]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

img = plt.imshow(model2, interpolation='none', cmap = cmap)

plt.show()

## iterate monte carlo moves
# for m in range(len(Temp)):
#     e0 = 0
#     e1 = 0
#     m0 = 0
#     m1 = 0
#     model = initial_state(N)
#     E = calc_energy(model)
#     M = calc_mag(model)
#
#     for i in range(equil_steps):
#         mc_moves(model, N, Temp[m])
#
#     for i in range(mc_steps):
#         mc_moves(model, N, Temp[m])
#         E = calc_energy(model)
#         M = calc_mag(model)
#
#
#         print(E,M,m,i)
#
#         # iterate
#         e0 = e0 + E
#         m0 = m0 + M
#         m1 = m1 + (M*M)
#         e1 = e1 + (E*E)
#
#         # change values in the arrays to the values
#         Energy[m]         = e0/(mc_steps*N*N)
#         Magnetization[m]  = m0/(mc_steps*N*N)
#         SpecificHeat[m]   = ( e1/mc_steps - e0*e0/(mc_steps*mc_steps) )/(N*Temp[m]*Temp[m])
#
#         Susceptibility[m] = ( m1/mc_steps - m0*m0/(mc_steps*mc_steps))





#######

# f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k');
#
# sp =  f.add_subplot(2, 2, 1 );
# plt.plot(Temp, Energy, 'o', color="red");
# plt.xlabel("Temperature (T)");
# plt.ylabel("Energy ", fontsize=20);
#
# sp =  f.add_subplot(2, 2, 2 );
# plt.plot(Temp, abs(Magnetization), 'o', color="blue");
# plt.xlabel("Temperature (T)");
# plt.ylabel("Magnetization ");
#
#
# sp =  f.add_subplot(2, 2, 3 );
# plt.plot(Temp, SpecificHeat, 'o', color="black");
# plt.xlabel("Temperature (T)");
# plt.ylabel("Specific Heat ");
#
#
# sp =  f.add_subplot(2, 2, 4 );
# plt.plot(Temp, Susceptibility, 'o', color="green");
# plt.xlabel("Temperature (T)");
# plt.ylabel("Susceptibility");
# plt.show()
# plt.legend(loc='best', fontsize=15);
