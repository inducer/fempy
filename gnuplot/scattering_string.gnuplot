# Data ------------------------------------------------------------------------
e=exp(1)
c2 = 1.2
c1 = 1

f(x)=e**(-x**2)
#f(x)=sin(x)/x

# Animation data --------------------------------------------------------------
t = 0
delta_t = 0.05
end_t = 10

# Problem definition ----------------------------------------------------------
shifted_f(x)=f(x+5)
# for x < 0
u1(x,t)=shifted_f(x-c1*t)
u2(x,t)=(c2-c1)/(c2+c1)*shifted_f(-c1*t-x)
# for x > 0
u3(x,t)=2*c2/(c2+c1)*shifted_f(c1/c2*(x-c2*t))
# on the whole line
u(x,t) = x < 0 ? u1(x,t)+u2(x,t) : u3(x,t)

# Setup -----------------------------------------------------------------------
max_y = abs( f(0.001) * 1.4 )
set samples 300
set yrange [-max_y:max_y]
set grid

plot u(x,t)
pause -1 "Press enter to start"

load "animator.gnuplot"
