import corner
import numpy as np
import matplotlib.pyplot as plt
import sys

#the directory where to get the data from
dir=sys.argv[1]

with open(dir+'/pew_output.dat') as file:
    data1 = [[digit for digit in line.split()] for line in file]


ndim=len(data1[1])-3
nsamples=len(data1)-1

names=data1[0][1:ndim]
print(names)

data=[]
for i in range(1,nsamples):
	data.append([])
	for j in range(0,ndim-1):
		data[i-1].append(float(data1[i][j]))


figure = corner.corner(np.array(data),labels=names,
                       quantiles=[0.16, 0.5, 0.84],smooth=1.0,bins=25,
                       show_titles=True, title_kwargs={"fontsize": 12},color='g')
                      
figure.savefig("cornerplot_"+dir+".pdf")

