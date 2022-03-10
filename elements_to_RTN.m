% Converts a 6x1 array of orbital elements [ a, ex, ey, i, R, u ]
% into RTN coordinates [ PR, PT, PN, VR, VT, VN ]
function [ dx_rtn ] = elements_to_RTN( oeC, oeD )
    
    % C ---> Chief or reference
    % D ---> Deputy satellite
    mu = 3.986004415e14;

    % Convert the orbital elements for SC1
    aC = oeC(1);
    wC = atan2( oeC(3), oeC(2) );
    eC = oeC(2)/cos(wC);
    iC = oeC(4);
    RC = oeC(5);
    uC = oeC(6);
    nuC = wrapToPi(uC - wC);

    % Convert the orbital elements for SC2
    aD = oeD(1);
    wD = atan2( oeD(3), oeD(2) );
    eD = oeD(2)/cos(wD);
    iD = oeD(4);
    RD = oeD(5);
    uD = oeD(6);
    nuD = wrapToPi(uD - wD);

    % Convert elements to ECI for Chief SC
    pC = aC*(1-eC^2);
    radC = pC / ( 1 + eC * cos(nuC) );
    pos_X = radC * cos( nuC );
    pos_Y = radC * sin( nuC );
    vel_X = -1 * sin(nuC) * sqrt( mu / pC );
    vel_Y = (eC + cos(nuC)) * sqrt( mu / pC );
    DCM_HN = rot313(RC,iC,wC);
    rC_eci = DCM_HN * [ pos_X pos_Y 0.0 ]' ;
    vC_eci = DCM_HN * [ vel_X vel_Y 0.0 ]' ;

    % Convert elements to ECI for Deputy SC
    pD = aD*(1-eD^2);
    radD = pD / ( 1 + eD * cos(nuD) );
    pos_X = radD * cos( nuD );
    pos_Y = radD * sin( nuD );
    vel_X = -1 * sin(nuD) * sqrt( mu / pD );
    vel_Y = (eD + cos(nuD)) * sqrt( mu / pD );
    DCM_HN = rot313(RD,iD,wD);
    rD_eci = DCM_HN * [ pos_X pos_Y 0.0 ]' ;
    vD_eci = DCM_HN * [ vel_X vel_Y 0.0 ]' ;
    
    % Angular velocity of RTN w.r.t. ECI, resolved in ECI frame % [rad/s]
    w_rtn_eci = cross( rC_eci, vC_eci) / norm(rC_eci)^2;
    
    % relative position and velocity resolved in ECI frame
    dr_eci = rD_eci-rC_eci;
    dv_eci = vD_eci-vC_eci-cross(w_rtn_eci,dr_eci);
    
    % R, T, and N basis vectors resolved in ECI frame
    R_eci = rC_eci / norm( rC_eci );
    N_eci = cross( rC_eci, vC_eci ) / norm(cross( rC_eci, vC_eci ));
    T_eci = cross( N_eci, R_eci );

    % Rotation matrix from ECI frame to RTN frame
    R_eci2rtn = [R_eci, T_eci, N_eci]';
    
    % Relative position and velocity resolved in RTN frame
    dr_rtn = R_eci2rtn*dr_eci;
    dv_rtn = R_eci2rtn*dv_eci;
    
    % Assembles RTN state of deputy
    dx_rtn = [dr_rtn;dv_rtn];
end