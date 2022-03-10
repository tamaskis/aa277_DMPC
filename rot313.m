function R313 = rot313(psi,theta,phi)
    sp = sin(psi);
    cp = cos(psi);
    st = sin(theta);
    ct = cos(theta);
    sf = sin(phi);
    cf = cos(phi);

    R313 = [ cf*cp-sf*ct*sp   cf*sp+sf*ct*cp   sf*st;
            -sf*cp-cf*ct*sp  -sf*sp+cf*ct*cp   cf*st;
             st*sp           -st*cp            ct]';
end