stacksize( 20000000 );
x = 0:0.1:1;
y = x;
z = ones(x'*y)-x'*y;
z = z($:-1:1,$:-1:1) .* z;
cmap = ( 0:0.01:1 )' * [ 1 1 1 ];
xset( "colormap", cmap );
Sgrayplot( x,y,z );
contour( x,y,z,10 );

