#!/usr/bin/env python3
import numpy as np
import mpi4py.MPI as MPI
import sys
import os
import matplotlib.pyplot as plt
if __name__ == "__main__":
    inputdirs = ['0_2','0.1_2','0.5_2','2_2','2_2LES','4_2']
    RFit = []
    XFit = []
    c=.16
    c1=2.63
    F = 2*np.pi**2*0.01*np.array([0,.1,0.5,2,2,4])
    symbols = ['o','^','s','d','o','x']
    for n,dir in enumerate(inputdirs):
        os.chdir('output_B'+dir+'level')
        print("Processing directory:", dir)
        data=np.load('results.npz'   )
        Z=data['CM'][:,2]
        X=(2*data['CM'][:,0])**2-1
        Y=(2*data['CM'][:,1])**2-1
        R=data['CM'][:,3]
        print(R[0])
        time = data['time']
        ZR= 1/(2*F[n])*np.array([1.5*R*np.log(R)+c1*(R-1)])
        if n>=2:
            plt.plot(R,Z,symbols[n],color='b', label='SOMAR')
            plt.plot(R, np.squeeze(ZR), '--', color='k', label='Turner 1957')
            plt.xlabel('R')
            plt.ylabel('z')
        
        RFit.append(np.polyfit(time, R**2-1,1))
        XFit.append(np.polyfit(Z,time,2))
        os.chdir('../')
    plt.savefig('Turner57_R_vs_Z.png')
    plt.close()
    plt.plot([f for f in F],[s[0]*np.pi for s in RFit],'o', color='b', label='SOMAR')
    plt.plot([f for f in F], [f for f in F],'--', color='k')
    
    plt.xlabel('Predicted total buoyancy')
    plt.ylabel('Calculated total buoyancy')
    plt.savefig('Turner57_R_vs_F.png')
    
    plt.close()
    
        


