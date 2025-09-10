import numpy as np


with open('TS-WAVE2D_FINERES/TS-WAVE2D_FINERES.geo.dat', 'r') as f:
    a = f.readlines()

conn = []
position = []

a = [x.strip() for x in a]
a = [x.strip('\n') for x in a]

# Read the connectivity matrix
for i,val in enumerate(a):
    #print(val)
    if val == 'ELEMENTS':
        print('IN ELEMENTS!')
        for ele in a[i+1:]:
            if ele != 'END_ELEMENTS':
                conn.append(ele)
            else:
                break
    if val == 'COORDINATES':
        print('IN COORDINATES!')
        for ele in a[i+1:]:
            if ele != 'END_COORDINATES':
                position.append(ele)
            else:
                break


conn = [x.split() for x in conn]
conn = [[int(x) for x in a[1:]] for a in conn]

position = [x.split() for x in position]
position = [[float(x) for x in a[1:]] for a in position]

conn = np.array(conn) # Connectivity matrix
position = np.array(position)

