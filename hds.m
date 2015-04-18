function [ gi ] = hds( wi )
T = 1;
fd = 5*10.^(-3);
wd = 2*pi*fd;
w0 = tan( wd*T/2 )
a1 = -( 2*(1-w0^2) )/( 1 + w0^2 + sqrt(2)*w0 );
a2 = ( 1 + w0^4 )/( 1 + w0^2 + sqrt(2)*w0 )^2;
c0 = ( 1 + a1 + a2 )/4;
gi = zeros( 1, length( wi ) );

for l = 3 : length( wi ) %saltiamo il transitorio.....
    gi( l ) = -a1*gi( l-1 ) -a2*gi( l-2 ) +c0*( wi( l ) +2*wi( l-1 )+wi( l-2 ) );
end
gi = gi( 3 : length( wi ) );
end

