f(x,y)=(1+x)*x*y**2*(1-x-y)**2/(x+y)/(x+(1-x-y))
g(x,y)=x+y>1?1/0:f(x,y)
splot [x=0:1][y=0:1] g(x,y)
