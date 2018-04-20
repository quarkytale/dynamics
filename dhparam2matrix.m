function T = dhparam2matrix(theta, d, a, alpha, i)
    r1 = [cos(theta(i)), -sin(theta(i))*cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i))];
    r2 = [sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i))];
    r3 = [0, sin(alpha(i)), cos(alpha(i)), d(i)];
    r4 = [0,0,0,1];
     
    T = [r1;r2;r3;r4];
end