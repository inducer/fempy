# source functions for the wave equation
# strauss p. 322

set terminal postscrip eps color
set isosamples 100
heaviside(x)=x>0?1:0
c=0.5
set hidden
set grid

# in 1d
set output "+wave_1d.eps"
splot [x=-5:5][t=-2:5] 1/(2*c)*heaviside(c**2*t**2-x**2)*sgn(t)
# in 2d
set output "+wave_2d.eps"
set zrange [-3000:3000]
splot [r=-0.01:0.01][t=-0.01:0.01] 1/(2*c)*heaviside(c**2*t**2-r**2)*1/sqrt(c**2*t**2-r**2)*sgn(t)
