function [ gi ] = hds( wi )
T = 1;
fd = 5*10.^(-3);
wd = 2*pi*fd;
w0 = tan( wd*T/2 );
a1 = -( 2*(1-w0^2) )/( 1 + w0^2 + sqrt(2)*w0 );
a2 = ( 1 + w0^4 )/( 1 + w0^2 + sqrt(2)*w0 )^2;
c0 = ( 1 + a1 + a2 )/4;
gi = zeros( 1, length( wi ) );
transient = 201;

%This part and the one Gian wrote (see commented part) are pretty much identical
a = [1 a1 a2];
b = [c0 2*c0 c0];
gi = filter(b, a, wi);
gi = gi(transient:length(gi));
% figure
% plot(1:length(gi), 10*log10(abs(gi)));


% for l = 3 : length( wi ) %saltiamo il transitorio.....
%     gi( l ) = -a1*gi( l-1 ) -a2*gi( l-2 ) +c0*( wi( l ) +2*wi( l-1 )+wi( l-2 ) );
% end
% gi = gi( transient : length( wi ) );
% figure
% plot(1:length(gi), 10*log10(abs(gi)));
% end

