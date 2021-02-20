#!/usr/bin/env python
# coding: utf-8

# ## SERIAL PLOTS
# 

# In[3]:



import numpy as np
from matplotlib import pyplot as plt


# In[39]:


# x is the plen * slen 
# y is the total time 

#VERSION 1
data_totaltime = np.array([
    [10000,0.0001 ],
    [100000,0.00046 ],
    [1000000, 0.00038],
    [10000000,0.0029 ],
    [1000000,0.00305],
    [100000000,0.02967 ],
    [1000000000,0.05085 ],
    [1000000000,0.28348],
    [10000000000,0.53645 ],
    [100000000000,5.23731 ],
    [100000000000,12.90995 ],
    [1000000000000, 138.59632]
])

data_matchtime = np.array([
    [10000,0.00003],
    [100000,0.0002],
    [1000000, 0.00019],
    [1000000,0.00194],
    [10000000,0.00177],
    [100000000,0.01949],
    [1000000000,0.04122],
    [1000000000,0.19065],
    [10000000000,0.44433],
    [100000000000, 4.33195],
    [100000000000,12.81745],
    [1000000000000, 137.68859]
])
    
x, y = data_totaltime.T

plt.figure(figsize=(25,12))
plt.plot(x,y, marker="o")
    
x, y = data_matchtime.T
plt.plot(x,y, marker="*")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.show()


# In[ ]:





# In[42]:


# x is the plen * slen 
# y is the total time 

#VERSION 2
data_totaltime = np.array([
    [10000,0.00010 ],
    [100000,0.00046],
    [1000000, 0.00038],
    [1000000,0.00296],
    [10000000,0.00287],
    [100000000, 0.02525],
    [1000000000,0.04979],
    [1000000000,0.24833],
    [10000000000,0.48612],
    [100000000000, 4.89723],
    [100000000000,13.19893],
    [1000000000000,145.08619]
])

data_matchtime = np.array([
    [10000, 0.00003],
    [100000,0.00019],
    [1000000, 0.00019],
    [1000000,0.00175],
    [10000000,0.00173],
    [100000000, 0.01547],
    [1000000000,0.03991],
    [1000000000,0.15494],
    [10000000000, 0.39396],
    [100000000000, 3.98295],
    [100000000000,13.10603],
    [1000000000000,144.16187]
])
x, y = data_totaltime.T
plt.figure(figsize=(25,12))
plt.plot(x,y, marker="o")

x, y = data_matchtime.T
plt.plot(x,y, marker="*")
plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.show()


# # Parallel code

# ### 4 CORES

# In[61]:


# Total times 4 cores
# Version 1 break 
# Version 1 task 
# Seriale 1

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial = np.array([
    [100000000, 0.00581 ],
    [1000000000,0.05085 ],
    [10000000000,0.53645 ],
    [100000000000,12.90995 ]

])

v1break = np.array([
    [100000000,0.00371],
    [1000000000,0.03421],
    [10000000000,0.32978],
    [100000000000,5.75996]
])

v1task= np.array([
    [100000000,0.01111],
    [1000000000,0.07921],
    [10000000000,0.75030],
    [100000000000,17.91615] # VALORE STRANO!!
])


x, y = serial.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial")

x, y = v1break.T
plt.plot(x,y, marker="*",  label="V1")

x, y = v1task.T
plt.plot(x,y, marker="+",  label="V1 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 8 cores V1")
plt.legend(loc="upper left")
plt.show()


# ### 8 CORES

# In[60]:


# Total times 8 cores
# Version 1 break 
# Version 1 task 
# Seriale 1

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial = np.array([
    [100000000, 0.00581 ],
    [1000000000,0.05085 ],
    [10000000000,0.53645 ],
    [100000000000,12.90995 ]

])

v1break = np.array([
    [100000000,0.01338],
    [1000000000,0.02391],
    [10000000000,0.23258],
    [100000000000,3.01167]
])

v1task= np.array([
    [100000000,0.01188],
    [1000000000,0.07670],
    [10000000000,0.73476],
    [100000000000,17.87497] #valore stranooooo
])


x, y = serial.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial")

x, y = v1break.T
plt.plot(x,y, marker="*",  label="V1")

x, y = v1task.T
plt.plot(x,y, marker="+",  label="V1 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 8 cores V1")
plt.legend(loc="upper left")
plt.show()


# ### 16 CORES

# In[57]:



# Total times 16 cores
# Version 1 break 
# Version 1 task 
# Seriale 1

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial = np.array([
   [100000000, 0.00581 ],
   [1000000000,0.05085 ],
   [10000000000,0.53645 ],
   [100000000000,12.90995 ]

])

v1break = np.array([
   [100000000, 0.00288],
   [1000000000,0.03225],
   [10000000000, 0.14792],
   [100000000000,1.45170]
])

v1task= np.array([
   [100000000, 0.00339],
   [1000000000,0.01545],
   [10000000000, 0.14558],
   [100000000000,1.44993]
])


x, y = serial.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial")

x, y = v1break.T
plt.plot(x,y, marker="*",  label="V1")

x, y = v1task.T
plt.plot(x,y, marker="+",  label="V1 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 16 cores V1")
plt.legend(loc="upper left")
plt.show()


# ### 24 CORES

# In[59]:


# Total times 24 cores

# Version 1 break 
# Version 1 task 
# Seriale 1

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial = np.array([
    [100000000, 0.00581 ],
    [1000000000,0.05085 ],
    [10000000000,0.53645 ],
    [100000000000,12.90995 ]

])

v1break = np.array([
    [100000000,0.00490],
    [1000000000,0.01441],
    [10000000000,0.13258],
    [100000000000,0.97405]
])

v1task= np.array([
    [100000000,0.00476],
    [1000000000,0.01405],
    [10000000000,0.12728],
    [100000000000,0.97600]
])


x, y = serial.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial")

x, y = v1break.T
plt.plot(x,y, marker="*",  label="V1")

x, y = v1task.T
plt.plot(x,y, marker="+",  label="V1 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 24 cores V1")
plt.legend(loc="upper left")
plt.show()


# ### 16 core V2

# In[58]:


#VERSIONE 2 PER AVERE UN CONFRONTO 

# Total times 16 cores
# Version 2 break 
# Version 2 task 
# Seriale 2

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial2 = np.array([
    [100000000, 0.00522],
    [1000000000,0.04979 ],
    [10000000000,0.48612 ],
    [100000000000,13.19893 ]

])

v2break = np.array([
    [100000000, 0.00619 ],
    [1000000000,0.03870 ],
    [10000000000, 0.38972 ],
    [100000000000,7.16041]
])

v2task= np.array([
    [100000000, 0.00210],
    [1000000000,0.06274],
    [10000000000, 0.58723],
    [100000000000,15.94546]
])


x, y = serial2.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial2")

x, y = v2break.T
plt.plot(x,y, marker="*",  label="V2")

x, y = v2task.T
plt.plot(x,y, marker="+",  label="V2 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 16 cores V2")
plt.legend(loc="upper left")
plt.show()


# #### 24 CORE V2

# In[62]:


#VERSIONE 2 PER AVERE UN CONFRONTO 

# Total times 16 cores
# Version 2 break 
# Version 2 task 
# Seriale 2

#Misurazioni: per S 10^5 P 10^3 - S 10^6 P 10^3 - S 10^7 P 10^3 - S 10^7 P 10^4


serial2 = np.array([
    [100000000, 0.00522],
    [1000000000,0.04979 ],
    [10000000000,0.48612 ],
    [100000000000,13.19893 ]

])

v2break = np.array([
    [100000000, 0.00269],
    [1000000000,0.04009 ],
    [10000000000, 0.46632],
    [100000000000,7.51481]
])

v2task= np.array([
    [100000000,0.00437 ],
    [1000000000,0.06373],
    [10000000000,0.59509],
    [100000000000,15.97497]
])


x, y = serial2.T
plt.figure(figsize=(25,8))
plt.plot(x,y, marker="o",  label="serial2")

x, y = v2break.T
plt.plot(x,y, marker="*",  label="V2")

x, y = v2task.T
plt.plot(x,y, marker="+",  label="V2 task")

plt.grid()
plt.xlabel("Slen x Plen")
plt.ylabel("Total Time [s]")
plt.title("Serial vs Parallel 24 cores V2")
plt.legend(loc="upper left")
plt.show()


# In[ ]:




