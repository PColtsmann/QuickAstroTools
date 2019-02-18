import numpy as np

"""

These functions rotate an orbit between a model frame and an external (Image) reference frame.
'rotate' functions move from external frame to model frame and 'derotate' functions move from model frame to external.

Those prefixed with 'exp_' have the explicit calculations written out, these are better suited for adjustments such that iteration along an axis can be performed.
Those without the 'exp_' prefix use numpy matrix multiplication.
Those with a 'g_' prefix are in the convention of the ZODIPIC code and Grant Kennedy's (drgmk) 'alma.alma.image' code for rotations. 

In all cases capital lettered variables denote the reference (Image) frame and lower cased variables denote model frame.

Either degrees or radians can be used for angles and either spherical polar (with elevation in place of polar angle) or cartesian coordinates can be used. In both cases the output will be in the same form as the input.

Arguments are as follows:

r/R     : radial distance of point of interest (either r/az/el can be used OR x/y/z to specify of point of interest)
az/AZ   : azimuthal angle of point of interest, defined East of North
el/EL   : elevation angle of point of interest
x/X     : X coordinate of point of interest
y/Y     : Y coordinate of point of interest
z/Z     : Z coordinate of point of interest
inc     : inclination of orbit
pos     : postion angle/longitude of ascending node of orbit, defined East of North
anom    : argument of periapsis/perifocus/pericentre of elliptical orbit
X0      : X coordinate offset of external frame relative to rotated model frame
Y0      : Y coordinate offset of external frame relative to rotated model frame
Z0      : Z coordinate offset of external frame relative to rotated model frame
DEG     : if left as 0 angular inputs are interpreted as in radians, otherwise (e.g. DEG = 1) inputs are interpreted as in degrees.

"""

def exp_rotate(R=0, AZ=0, EL=0, X=0, Y=0, Z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):

    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        AZ = np.deg2rad(AZ)
        EL = np.deg2rad(EL)     
    if R or AZ or EL:
        THE = np.pi/2 -  EL
        X = R*np.cos(AZ)*np.sin(THE) - X0
        Y = R*np.sin(AZ)*np.sin(THE) - Y0
        Z = R*np.cos(THE)             - Z0
    elif X or Y or Z:
        X = X - X0
        Y = Y - Y0
        Z = Z - Z0
        
    c0 = np.cos(pos)
    s0 = np.sin(pos)
    c1 = np.cos(inc)
    s1 = np.sin(inc)
    c2 = np.cos(anom)
    s2 = np.sin(anom)
    
    x = (-s0*c1*s2 + c0*c2)*X + (c0*c1*s2 + s0*c2)*Y + (s1*s2)*Z
    y = (-s0*c1*c2 - c0*s2)*X + (c0*c1*c2 - s0*s2)*Y + (s1*c2)*Z
    z =             (s0*s1)*X +           (-c0*s1)*Y +    (c1)*Z
    
    rxy2 = x**2 + y**2
    rxy = np.sqrt(rxy2)
    r = np.sqrt(rxy2 + z**2)
    az = np.arctan2(y,x)
    el = np.arctan2(z,rxy)
    
    if DEG:
        az = np.rad2deg(az)
        el = np.rad2deg(el)
    if R or AZ or EL:
        return [r,az,el]          
    elif X or Y or Z:
        return [x,y,z]
    
def rotate(R=0, AZ=0, EL=0, X=0, Y=0, Z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        AZ = np.deg2rad(AZ)
        EL = np.deg2rad(EL)     
    if R or AZ or EL:
        THE = np.pi/2 -  EL
        X = R*np.cos(AZ)*np.sin(THE) - X0
        Y = R*np.sin(AZ)*np.sin(THE) - Y0
        Z = R*np.cos(THE)             - Z0
    elif X or Y or Z:
        X = X - X0
        Y = Y - Y0
        Z = Z - Z0
        
    IMPLANE = np.array([[X],[Y],[Z]])

    P3pos = np.array([[np.cos(pos),-np.sin(pos),0],
                      [np.sin(pos), np.cos(pos),0],
                      [      0     ,      0    ,1]])
    P2inc = np.array([[1,    0     ,      0      ],
                      [0,np.cos(inc),-np.sin(inc)],
                      [0,np.sin(inc), np.cos(inc)]])
    P1ano = np.array([ [np.cos(anom),-np.sin(anom),0],
                       [np.sin(anom), np.cos(anom),0],
                       [     0      ,      0      ,1]])
    P321 = np.matmul(P3pos,np.matmul(P2inc,P1ano))
    MODPLANE = np.matmul(np.transpose(P321),IMPLANE) 

    x = MODPLANE[0][0]
    y = MODPLANE[1][0]
    z = MODPLANE[2][0]
    
    rxy2 = x**2 + y**2
    rxy = np.sqrt(rxy2)
    r = np.sqrt(rxy2 + z**2)
    az = np.arctan2(y,x)
    el = np.arctan2(z,rxy)
    
    if DEG:
        az = np.rad2deg(az)
        el = np.rad2deg(el)
    if R or AZ or EL:
        return [r,az,el]          
    elif X or Y or Z:
        return [x,y,z]
  
def exp_derotate(r=0, az=0, el=0, x=0, y=0, z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        az = np.deg2rad(az)
        el = np.deg2rad(el)     
    if r or az or el:
        the = np.pi/2 -  el
        x = r*np.cos(az)*np.sin(the)
        y = r*np.sin(az)*np.sin(the)
        z = r*np.cos(the)            
    elif x or y or z:
        x = x
        y = y
        z = z  

    c0 = np.cos(pos)
    s0 = np.sin(pos)
    c1 = np.cos(inc)
    s1 = np.sin(inc)
    c2 = np.cos(anom)
    s2 = np.sin(anom)
    
    X = (c0*c2 - c1*s0*s2)*x + (-c0*s2 - c1*s0*c2)*y +  (s1*s0)*z + X0
    Y = (s0*c2 + c1*c0*s2)*x +  (c1*c0*c2 - s0*s2)*y + (-s1*c0)*z + Y0
    Z =            (s1*s2)*x +             (s1*c2)*y +     (c1)*z + Z0

    RXY2 = X**2 + Y**2
    RXY = np.sqrt(RXY2)
    R = np.sqrt(RXY2 + Z**2)
    AZ = np.arctan2(Y,X)
    EL = np.arctan2(Z,RXY)
    
    if DEG:
        AZ = np.rad2deg(AZ)
        EL = np.rad2deg(EL)
    if r or az or el:
        return [R,AZ,EL]          
    elif x or y or z:
        return [X,Y,Z]
    

def derotate(r=0, az=0, el=0, x=0, y=0, z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        az = np.deg2rad(az)
        el = np.deg2rad(el)     
    if r or az or el:
        the = np.pi/2 -  el
        x = r*np.cos(az)*np.sin(the)
        y = r*np.sin(az)*np.sin(the)
        z = r*np.cos(the)             
    elif x or y or z:
        x = x
        y = y
        z = z
        
    IMPLANE = np.array([[x],[y],[z]])

    P3pos = np.array([ [np.cos(pos),-np.sin(pos),0],
                       [np.sin(pos), np.cos(pos),0],
                       [     0     ,       0    ,1]])
    P2inc = np.array([[1,    0      ,      0     ],
                      [0,np.cos(inc),-np.sin(inc)],
                      [0,np.sin(inc), np.cos(inc)]])
    P1ano = np.array([ [np.cos(anom),-np.sin(anom),0],
                       [np.sin(anom), np.cos(anom),0],
                       [     0      ,       0     ,1]])
    P321 = np.matmul(P3pos,np.matmul(P2inc,P1ano))
    MODPLANE = np.matmul(P321,IMPLANE) 
    
    X = MODPLANE[0][0] + X0
    Y = MODPLANE[1][0] + Y0
    Z = MODPLANE[2][0] + Z0
    
    RXY2 = X**2 + Y**2
    RXY = np.sqrt(RXY2)
    R = np.sqrt(RXY2 + Z**2)
    AZ = np.arctan2(Y,X)
    EL = np.arctan2(Z,RXY)
    
    if DEG:
        AZ = np.rad2deg(AZ)
        EL = np.rad2deg(EL)
    if r or az or el:
        return [R,AZ,EL]          
    elif x or y or z:
        return [X,Y,Z]
        
"""
The following functions use the convention of the ZODIPIC code and Grant Kennedy's (drgmk) 'alma.alma.image' code for rotations.
"""

def g_exp_rotate(R=0, AZ=0, EL=0, X=0, Y=0, Z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    """Image to Model"""
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        AZ = np.deg2rad(AZ)
        EL = np.deg2rad(EL)     
    if R or AZ or EL:
        THE = np.pi/2 -  EL
        X = R*np.cos(AZ)*np.sin(THE) - X0
        Y = R*np.sin(AZ)*np.sin(THE) - Y0
        Z = R*np.cos(THE)             - Z0
    elif X or Y or Z:
        X = X - X0
        Y = Y - Y0
        Z = Z - Z0
        
    c0 = np.cos(-pos)
    s0 = np.sin(-pos)
    c1 = np.cos(-inc)
    s1 = np.sin(-inc)
    c2 = np.cos(-anom-np.pi/2)
    s2 = np.sin(-anom-np.pi/2)
    
    x = (-s0*c1*s2 + c0*c2)*Y + (c0*c1*s2 + s0*c2)*X + (s1*s2)*Z
    y = (-s0*c1*c2 - c0*s2)*Y + (c0*c1*c2 - s0*s2)*X + (s1*c2)*Z
    z =             (s0*s1)*Y +           (-c0*s1)*X +    (c1)*Z
    
    rxy2 = x**2 + y**2
    rxy = np.sqrt(rxy2)
    r = np.sqrt(rxy2 + z**2)
    az = np.arctan2(y,x)
    el = np.arctan2(z,rxy)
    
    if DEG:
        az = np.rad2deg(az)
        el = np.rad2deg(el)
    if R or AZ or EL:
        return [r,az,el]          
    elif X or Y or Z:
        return [x,y,z]
    
def g_rotate(R=0, AZ=0, EL=0, X=0, Y=0, Z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        AZ = np.deg2rad(AZ)
        EL = np.deg2rad(EL)     
    if R or AZ or EL:
        THE = np.pi/2 -  EL
        X = R*np.cos(AZ)*np.sin(THE) - X0
        Y = R*np.sin(AZ)*np.sin(THE) - Y0
        Z = R*np.cos(THE)            - Z0
    elif X or Y or Z:
        X = X - X0
        Y = Y - Y0
        Z = Z - Z0
    IMPLANE = np.array([[Y],[X],[Z]])

    P3pos = np.array([ [np.cos(-pos),-np.sin(-pos),0],
                       [np.sin(-pos), np.cos(-pos),0],
                       [      0     ,      0      ,1]])
    P2inc = np.array([[1,    0       ,      0      ],
                      [0,np.cos(-inc),-np.sin(-inc)],
                      [0,np.sin(-inc), np.cos(-inc)]])
    P1ano = np.array([ [np.cos(-anom-np.pi/2),-np.sin(-anom-np.pi/2),0],
                       [np.sin(-anom-np.pi/2), np.cos(-anom-np.pi/2),0],
                       [     0               ,      0               ,1]])
    P321 = np.matmul(P3pos,np.matmul(P2inc,P1ano))
    MODPLANE = np.matmul(np.transpose(P321),IMPLANE) 

    x = MODPLANE[0][0]
    y = MODPLANE[1][0]
    z = MODPLANE[2][0]
    
    rxy2 = x**2 + y**2
    rxy = np.sqrt(rxy2)
    r = np.sqrt(rxy2 + z**2)
    az = np.arctan2(y,x)
    el = np.arctan2(z,rxy)
    
    if DEG:
        az = np.rad2deg(az)
        el = np.rad2deg(el)
    if R or AZ or EL:
        return [r,az,el]          
    elif X or Y or Z:
        return [x,y,z]
  
def g_exp_derotate(r=0, az=0, el=0, x=0, y=0, z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):

    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        az = np.deg2rad(az)
        el = np.deg2rad(el)     
    if r or az or el:
        the = np.pi/2 -  el
        x = r*np.cos(az)*np.sin(the)
        y = r*np.sin(az)*np.sin(the)
        z = r*np.cos(the)            
    elif x or y or z:
        x = x
        y = y
        z = z  

    c0 = np.cos(-pos)
    s0 = np.sin(-pos)
    c1 = np.cos(-inc)
    s1 = np.sin(-inc)
    c2 = np.cos(-anom-np.pi/2)
    s2 = np.sin(-anom-np.pi/2)
    
    Y = (c0*c2 - c1*s0*s2)*x + (-c0*s2 - c1*s0*c2)*y +  (s1*s0)*z + Y0
    X = (s0*c2 + c1*c0*s2)*x +  (c1*c0*c2 - s0*s2)*y + (-s1*c0)*z + X0
    Z =            (s1*s2)*x +             (s1*c2)*y +     (c1)*z + Z0

    RXY2 = X**2 + Y**2
    RXY = np.sqrt(RXY2)
    R = np.sqrt(RXY2 + Z**2)
    AZ = np.arctan2(Y,X)
    EL = np.arctan2(Z,RXY)
    
    if DEG:
        AZ = np.rad2deg(AZ)
        EL = np.rad2deg(EL)
    if r or az or el:
        return [R,AZ,EL]          
    elif x or y or z:
        return [X,Y,Z]
    
def g_derotate(r=0, az=0, el=0, x=0, y=0, z=0, inc=0, pos=0, anom=0, X0=0, Y0=0, Z0=0, DEG=0):
    
    if DEG:
        inc = np.deg2rad(inc)
        pos = np.deg2rad(pos)
        anom= np.deg2rad(anom)
        az = np.deg2rad(az)
        el = np.deg2rad(el)     
    if r or az or el:
        the = np.pi/2 -  el
        x = r*np.cos(az)*np.sin(the)
        y = r*np.sin(az)*np.sin(the)
        z = r*np.cos(the)             
    elif x or y or z:
        x = x
        y = y
        z = z
        
    IMPLANE = np.array([[x],[y],[z]])

    P3pos = np.array([ [np.cos(-pos),-np.sin(-pos),0],
                       [np.sin(-pos), np.cos(-pos),0],
                       [      0     ,      0      ,1]])
    P2inc = np.array([[1,    0       ,      0      ],
                      [0,np.cos(-inc),-np.sin(-inc)],
                      [0,np.sin(-inc), np.cos(-inc)]])
    P1ano = np.array([ [np.cos(-anom-np.pi/2),-np.sin(-anom-np.pi/2),0],
                       [np.sin(-anom-np.pi/2), np.cos(-anom-np.pi/2),0],
                       [     0               ,      0               ,1]])
    P321 = np.matmul(P3pos,np.matmul(P2inc,P1ano))
    MODPLANE = np.matmul(P321,IMPLANE) 
    
    Y = MODPLANE[0][0] + Y0
    X = MODPLANE[1][0] + X0
    Z = MODPLANE[2][0] + Z0
    
    RXY2 = X**2 + Y**2
    RXY = np.sqrt(RXY2)
    R = np.sqrt(RXY2 + Z**2)
    AZ = np.arctan2(Y,X)
    EL = np.arctan2(Z,RXY)
    
    if DEG:
        AZ = np.rad2deg(AZ)
        EL = np.rad2deg(EL)
    if r or az or el:
        return [R,AZ,EL]          
    elif x or y or z:
        return [X,Y,Z]
