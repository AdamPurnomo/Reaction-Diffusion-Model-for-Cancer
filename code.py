import numpy as np
import matplotlib.pyplot as plt

rootA = "Results/A/"
rootB = "Results/B/"
rootC = "Results/C/"
rootD = "Results/C/"

x,y = np.mgrid[-2.5:2.5:0.1, -2.5:2.5:0.1]
length, width = x.shape

H = np.zeros((100000,width,length))
C = np.zeros((100000,width,length))
L = np.zeros((100000,width,length))

#initialize cancer cell, healthy cell and acid concentration
H[0] = np.random.uniform(low=4.91e7, high=4.99e7, size=(width,length))
C[0] = 3e7*np.exp(-5*np.sqrt(x**2 + y**2))

# Reaction Parameters
r_H = 1e-6
r_C = 1e-6
r_L = 2.2e-17
d_H = 2*r_H
d_C = 3*r_C

k_H = 5e7
k_C = 5e7
alpha = 0.0
beta = 0.5
gamma = 1.1e-4

#Diffusion Parameters
D_C = 2e-10
D_L = 5e-6 
T = 0


#time and space parameters
delta_space = 0.1
delta_time = 300


#time init
t = 1

def calculate(t):
    #reaction terms
    reac_H = r_H*H[t]*(1 - (H[t]/k_H) - (alpha*C[t]/k_C)) - beta*L[t]*H[t] - d_H*T[t]*H[t]*(1 - (H[t]/k_H)) 
    reac_C = r_C*C[t]*(1 - (C[t]/k_C) - (alpha*H[t]/k_H)) - d_C*T[t]*C[t]*(1 - (C[t]/k_C))
    reac_L = r_L*C[t] - gamma*L[t]

    ##Diffusion terms

    #healthy cell
    diff_H = 0

    #cancer cell
    kappa = D_C*(1 - H[t]/k_H) #diffusion coefficient of cancer as a function of healthy cell
    diff_C = np.zeros((width,length))
    for i in range(width):
        for j in range(length):
            if((i==0) and (j==0)):
                diff_C[i,j] = (1/(delta_space**2))*((kappa[i,j] + kappa[i,j+1])*(C[t,i,j+1] - C[t,i,j]) + ((kappa[i,j] + kappa[i+1,j])*(C[t,i+1,j] - C[t,i,j])))
            elif((i==0) and (j==length-1)):
                diff_C[i,j] = (1/(delta_space**2))*((kappa[i,j] + kappa[i,j-1])*(C[t,i,j-1] - C[t,i,j]) + ((kappa[i,j] + kappa[i+1,j])*(C[t,i+1,j] - C[t,i,j])))
            elif((i==width-1) and (j==length-1)):
                diff_C[i,j] = (1/(delta_space**2))*((kappa[i,j] + kappa[i,j-1])*(C[t,i,j-1] - C[t,i,j]) + ((kappa[i,j] + kappa[i-1,j])*(C[t,i-1,j] - C[t,i,j])))
            elif((i==width-1) and (j==0)):
                diff_C[i,j] = (1/(delta_space**2))*((kappa[i,j] + kappa[i,j+1])*(C[t,i,j+1] - C[t,i,j]) + ((kappa[i,j] + kappa[i-1,j])*(C[t,i-1,j] - C[t,i,j])))
            elif(j==0):
                diff_C[i,j] = (1/(2*delta_space**2))*(2*(kappa[i,j+1] + kappa[i,j])*(C[t,i,j+1] - C[t,i,j]) + (kappa[i+1,j] + kappa[i,j])*(C[t,i+1,j] - C[t,i,j]) - (kappa[i,j] + kappa[i-1,j])*(C[t,i,j] - C[t,i-1,j]))
            elif(j==length-1):
                diff_C[i,j] = (1/(2*delta_space**2))*(2*(kappa[i,j] + kappa[i,j-1])*(C[t,i,j-1] - C[t,i,j]) + (kappa[i+1,j] + kappa[i,j])*(C[t,i+1,j] - C[t,i,j]) - (kappa[i,j] + kappa[i-1,j])*(C[t,i,j] - C[t,i-1,j]))
            elif(i==0):
                diff_C[i,j] = (1/(2*delta_space**2))*((kappa[i,j+1] + kappa[i,j])*(C[t,i,j+1] - C[t,i,j]) - (kappa[i,j] + kappa[i,j-1])*(C[t,i,j] - C[t,i,j-1]) + 2*(kappa[i,j] + kappa[i+1,j])*(C[t,i+1,j] - C[t,i,j]))
            elif(i==width-1):
                diff_C[i,j] = (1/(2*delta_space**2))*((kappa[i,j+1] + kappa[i,j])*(C[t,i,j+1] - C[t,i,j]) - (kappa[i,j] + kappa[i,j-1])*(C[t,i,j] - C[t,i,j-1]) + 2*(kappa[i-1,j] + kappa[i,j])*(C[t,i-1,j] - C[t,i,j]))
            else:
                diff_C[i,j] = (1/(2*delta_space**2))*((kappa[i,j+1] + kappa[i,j])*(C[t,i,j+1] - C[t,i,j]) - (kappa[i,j] + kappa[i,j-1])*(C[t,i,j] - C[t,i,j-1]) + (kappa[i+1,j] + kappa[i,j])*(C[t,i+1,j] - C[t,i,j]) - (kappa[i,j] + kappa[i-1,j])*(C[t,i,j] - C[t,i-1,j]))




            

    #acid concetration
    diff_L = np.zeros((width,length))
    #only iterating outside boundary region
    for i in range(width):
        for j in range(length):
            if((i==0) and (j==0)):
                diff_L[i,j] = (2/(delta_space**2))*(L[t,i,j+1] - L[t,i,j] + L[t,i+1,] - L[t,i,j])
            elif((i==0) and (j==length-1)):
                diff_L[i,j] = (2/(delta_space**2))*(L[t,i,j-1] - L[t,i,j] + L[t,i+1,] - L[t,i,j])
            elif((i==width-1) and (j==length-1)):
                diff_L[i,j] = (2/(delta_space**2))*(L[t,i,j-1] - L[t,i,j] + L[t,i-1,] - L[t,i,j])
            elif((i==width-1) and (j==0)):
                diff_L[i,j] = (2/(delta_space**2))*(L[t,i,j+1] - L[t,i,j] + L[t,i-1,] - L[t,i,j])
            elif(j==0):
                diff_L[i,j] = (1/(delta_space**2))*(2*(L[t,i,j+1] - L[t,i,j]) + (L[t,i+1,j] - 2*L[t,i,j] + L[t,i-1,j]))
            elif(j==length-1):
                diff_L[i,j] = (1/(delta_space**2))*(2*(L[t,i,j-1] - L[t,i,j]) + (L[t,i+1,j] - 2*L[t,i,j] + L[t,i-1,j]))
            elif(i==0):
                diff_L[i,j] = (1/(delta_space**2))*((L[t,i,j+1] - 2*L[t,i,j] + L[t,i,j-1]) + 2*(L[t,i+1,j] - L[t,i,j]))
            elif(i==width-1):
                diff_L[i,j] = (1/(delta_space**2))*((L[t,i,j+1] - 2*L[t,i,j] + L[t,i,j-1]) + 2*(L[t,i-1,j] - L[t,i,j]))
            else:
                diff_L[i,j] = (1/(delta_space**2))*((L[t,i,j+1] - 2*L[t,i,j] + L[t,i,j-1]) + (L[t,i+1,j] - 2*L[t,i,j] + L[t,i-1,j]))


    C[t+1] = C[t] + delta_time*(D_C*diff_C + reac_C)
    H[t+1] = H[t] + delta_time*(diff_H + reac_H)
    L[t+1] = L[t] + delta_time*(D_L*diff_L + reac_L)
        
    #make sure that it does not go negative
    invalid_idx = np.where(C[t+1]<0)
    C[t+1,invalid_idx] = 0

from mpl_toolkits.axes_grid1 import make_axes_locatable

def visualize(time, treament):
    plt.clf()
    f,ax = plt.subplots(1,3)
    
    im1 = ax[0].imshow(C[time], cmap="Reds", vmin=0, vmax=5e7)
    ax[0].set_title("Cancerous Cells", pad=15)
    divider = make_axes_locatable(ax[0])
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    cb1 = f.colorbar(im1,cax=cax1)
    cb1.ax.yaxis.set_offset_position('left')
    cb1.update_ticks()
    
    im2 = ax[1].imshow(L[time], cmap="Blues", vmin=0, vmax=9e-6)
    ax[1].set_title("Acid Concentration", pad=15)
    divider = make_axes_locatable(ax[1])
    cax2 = divider.append_axes("right", size="5%", pad=0.05)
    cb2 = f.colorbar(im2,cax=cax2)
    cb2.ax.yaxis.set_offset_position('left')
    cb2.update_ticks()
    
    im3 = ax[2].imshow(H[time], cmap="Greens", vmin=0, vmax=5e7)
    ax[2].set_title("Normal Cells", pad=15)
    divider = make_axes_locatable(ax[2])
    cax3 = divider.append_axes("right", size="5%", pad=0.05)
    cb3 = f.colorbar(im3,cax=cax3)
    cb3.ax.yaxis.set_offset_position('left')
    cb3.update_ticks()
    
    if(treament=="A"):
        root = rootA
    elif(treament=="B"):
        root = rootB
    elif(treament=="C"):
        root = rootC
    elif(treament=="D"):
        root = rootD

    day = time//288
    f.suptitle('Day ' + str(day), y=0.8)
    f.tight_layout()
    f.savefig(root + "Day" + str(day) + ".png",bbox_inches='tight')
    
    plt.show()

chemo = True
surgery = True
tday = 288 #num of iteration for one day period
if(chemo==True):
    T = np.zeros((100000))

    t_start = 30*tday
    t_on1 = 44*tday
    t_off1 = 58*tday
    t_on2 = 72*tday
    t_off2 = 86*tday
    t_on3 = 72*tday
    t_off3 = 86*tday
    t_on4 = 100*tday
    t_off4 = 114*tday
    t_on5 = 128*tday
    t_off5 = 142*tday
    t_on6 = 156*tday

    T[t_start:t_on1] = 1
    T[t_off1:t_on2] = 1
    T[t_off2:t_on3] = 1
    T[t_off3:t_on4] = 1
    T[t_off4:t_on5] = 1
    T[t_off5:t_on6] = 1

## Main loop
day = 0
for t in range(100000):
    if(t%288==0):
        day+=1
        visualize(t, "D")
    if(surgery==True):
        t_surgery = 30*tday
        if(t==t_surgery):
            C[t] = (0.01)*C[t]
            L[t] = np.zeros((width,length))
    if(day==301):
        break
    calculate(t)
    