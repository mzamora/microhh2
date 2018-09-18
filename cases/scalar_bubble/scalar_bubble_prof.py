import numpy

# set geometry
zsize=1.
kmax=64

# define variables
z=numpy.zeros(kmax)
u=numpy.zeros(kmax)
myscalar=numpy.zeros(kmax)

for k in range(kmax):
  z[k]=k*zsize/(kmax-1)
  u[k]=0.1

#write data
profile=open('scalar_bubble.prof','w')
profile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z','u','myscalar'))
for k in range(kmax):
  profile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(z[k], u[k],  myscalar[k]))
profile.close()
