plot u(x,t),u1(x,t),u2(x,t)
! python -c "import time;time.sleep(0.05)"
t = t + delta_t
if ( t < end_t ) reread
