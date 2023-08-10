# ICBM (Intercontinental Ballistic Missile) trajectory simulation
# www.overfitting.net
# https://www.overfitting.net/2023/08/oppenheimer-i-am-become-death-destroyer.html


library(ggmap)  # map_data()
library(data.table)


# BITMAP DRAWING FUNCTIONS

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}

DrawEllip = function(img, x0, y0, a, b, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja elipse de centro (x0,y0) y radios a y b
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    # Aquí no redondeamos para tener más precisión en la división
    if (fill) {
        indices=which( ((row(img)-x0)/a)^2 + ((col(img)-y0)/b)^2 < 1 )
    } else {
        indices=which( ((row(img)-x0)/(a+thick/2))^2 + ((col(img)-y0)/(b+thick/2))^2 <  1 &
                       ((row(img)-x0)/(a-thick/2))^2 + ((col(img)-y0)/(b-thick/2))^2 >= 1 )
    }
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

DrawCircle = function(img, x0, y0, r, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja círculo de centro (x0,y0) y radio r
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    img=DrawEllip(img, x0, y0, r, r, inc, val, fill, thick)
    
    return(img)
}

DrawGradCircle = function(img, x0, y0, r, valmin=0, valmax=1, gamma=3.3) {
    indices=which( ((row(img)-x0)/r)^2 + ((col(img)-y0)/r)^2 < 1 )
    imgtmp=(valmax-valmin) * ((r-((row(img)-x0)^2 + (col(img)-y0)^2)^0.5)/r)^(1/gamma)
    img[indices]=imgtmp[indices]+valmin

    return(img)
}

# Por Carlos Gil Bellosta
indices.drawline = function(x0, y0, x1, y1) {
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (y0 == y1) return(cbind(x0:x1, y0)) # Recta de m=0 o un punto
    if (abs(x1 - x0) >= abs(y1 - y0)) { # Recta de 0 < |m| <= 1
        m = (y1 - y0) / (x1 - x0)
        cbind(x0:x1, round(y0 + m * ((x0:x1) - x0)))
    } else indices.drawline(y0, x0, y1, x1)[, 2:1]  # Recta de |m| > 1
    # Llamada traspuesta recursiva y traspuesta
}

DrawLine = function(img, x0, y0, x1, y1, inc=TRUE, val=1) {
    # Dibuja recta desde (x0,y0)-(x1,y1)
    # Por defecto método no destructivo y con valor=1
    indices=indices.drawline(x0, y0, x1, y1)
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

DrawPoint = function(img, x0, y0, inc=TRUE, val=1) {
    # Dibuja punto en (x0,y0)
    # Por defecto método no destructivo y con valor=1
    img=DrawLine(img, x0, y0, x0, y0, inc, val)
    
    return(img)
}

ShowBitmap = function(img, trunc=TRUE, gamma=1, chan=2, interpolate=F) {
    # Muestra bitmap en pantalla
    # Solo si trunc=F y la imagen excede de 1 se reescala a 1
    # Si no es monocromo se muestra el canal chan (por defecto G)
    if (length(dim(img))>2) img=img[,,chan]
    img[img<0]=0
    if (trunc) img[img>1]=1
    plot(as.raster(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), max=1),
         interpolate=interpolate)
}

SaveBitmap = function(img, name, trunc=TRUE, gamma=1) {
    # Guarda bitmap en formato PNG
    # Solo si trunc=FALSE y la imagen excede de 1 se reescala a 1
    library(png)
    img[img<0]=0
    if (trunc) img[img>1]=1
    if (tolower(substr(name, nchar(name)-3, nchar(name))) != ".png") name=paste0(name,".png")
    writePNG(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name)
}

LoadBitmap = function(name, chan=2) {
    # Lee bitmap en formato PNG
    # Si no es monocromo se carga el canal chan (por defecto G)
    library(png)
    img=readPNG(name)
    if (length(dim(img))>2) img=img[,,chan]
    
    return(t(img[nrow(img):1,]))
}


# COORDINATES CONVERSION

# Polar to XYZ
polar2x = function(r, phi, theta) { r*cos(theta)*sin(phi)}
polar2y = function(r, phi, theta) { r*sin(theta)}
polar2z = function(r, phi, theta) {-r*cos(theta)*cos(phi)}

# XYZ to Polar (a bit inefficient)
xyz2r     = function(x, y, z) {(x*x+y*y+z*z)^0.5}
xyz2phi   = function(x, y, z) {
    # if (x==0) {if (z>0) return (pi) else return (0)}
    acos( -z/( (x*x+y*y+z*z)^0.5 * cos(asin(y/(x*x+y*y+z*z)^0.5)) ) ) * sign(x)
    }
xyz2theta = function(x, y, z) {asin(y/(x*x+y*y+z*z)^0.5)}




# These functions expect a dataframe with 3 columns (x,y,z)
translateXYZ = function(df, dx=0, dy=0, dz=0) {  # translation
    df[,1]=df[,1]+dx
    df[,2]=df[,2]+dy
    df[,3]=df[,3]+dz
    return(df)
}

scaleXYZ = function(df, sx=1, sy=1, sz=1) {  # scale
    df[,1]=df[,1]*sx
    df[,2]=df[,2]*sy
    df[,3]=df[,3]*sz
    return(df)
}

rotateX = function(df, theta=0) {  # rotation around X axis
    df$xr = df[,1]
    df$yr = cos(theta)*df[,2] - sin(theta)*df[,3]
    df$zr = sin(theta)*df[,2] + cos(theta)*df[,3]
    df=df[,4:6]  # keep only rotated coords
    colnames(df)=c('x','y','z')
    return(df)
}

rotateY = function(df, theta=0) {  # rotation around Y axis
    df$xr = cos(theta)*df[,1] + sin(theta)*df[,3]
    df$yr = df[,2]    
    df$zr =-sin(theta)*df[,1] + cos(theta)*df[,3]
    df=df[,4:6]  # keep only rotated coords
    colnames(df)=c('x','y','z')
    return(df)
}

rotateZ = function(df, theta=0) {  # rotation around Z axis
    df$xr = cos(theta)*df[,1] - sin(theta)*df[,2]
    df$yr = sin(theta)*df[,1] + cos(theta)*df[,2]
    df$zr = df[,3]
    df=df[,4:6]  # keep only rotated coords
    colnames(df)=c('x','y','z')
    return(df)
}


PlotTraj = function() {  # lazy function using global variables
    trajplottmp$factor=f/trajplottmp$z
    trajplottmp$xp=round(trajplottmp$x*trajplottmp$factor + NCOLDIV2)
    trajplottmp$yp=round(trajplottmp$y*trajplottmp$factor + NROWDIV2)
    trajplottmp=trajplottmp[trajplottmp$xp>=1 & trajplottmp$xp<=DIMX &
                            trajplottmp$yp>=1 & trajplottmp$yp<=DIMY]
    # (xp, yp) is a 3D to 2D projection rounded and clipped
    ifboom=trajplottmp[trajplottmp$boom==1]
    if (nrow(ifboom)==1) {  # nuke to be plotted
        scale=max(3-abs(frame-Npoints[i]-1), 1)  # scale? 1-2-3-2-1-1-1-1...
        img=DrawCircle(img, ifboom$xp, ifboom$yp, RADIUSNUKE*scale,
                       inc=FALSE, fill=TRUE, val=GRAYNUKE)
    }
    img[cbind(trajplottmp$xp, trajplottmp$yp)]=GRAYICBM*trajplottmp$grayscale
    
    return (img)
}


################################################################################

# Clear all warnings()
assign("last.warning", NULL, envir = baseenv())


# DRAW EARTH

# PHYSICAL PARAMETERS
Rearth=6371.23  # Earth average radius (km)
Rmoon=1737.4  # Moon average radius (km)
dearthmoon=385000  # centre to centre Earth to Moon distance (km)
# dobserver.iss=408  # ISS average altitude (km)
# dobserver.moon=dearthmoon-Rearth-Rmoon  # observation point to Earth surface distance (km)
dz=dearthmoon-Rmoon  # observation point to Earth surface distance (km)
thetamax=acos(Rearth/dz)
distmax=(dz^2 - Rearth^2)^0.5  # max distance to visible points


# ANIMATION PARAMETERS
DIMX=1920  # Full HD animation
DIMY=1080  # 1920 x 1080 pixels
NCOLDIV2=round(DIMX/2)
NROWDIV2=round(DIMY/2)
TH=0.7  # border between Earth globe and image limits
RADIUS=min(NCOLDIV2, NROWDIV2)*TH
MAXH=3000  # max ICBM altitude from longest distance ICBM (km)
RADIUSNUKE=4  # boom circle mark
GRAYGLOBEMAX=0.35
GRAYGLOBEMIN=0.05
GRAYMAP=0.6
GRAYICBM=1
GRAYNUKE=1

# Calculate focal length to fit Earth in the final image and FOV
f=min(NCOLDIV2,NROWDIV2)*TH*
    (dz-Rearth*cos(thetamax))/(Rearth*sin(thetamax))
FOV=(pi-2*thetamax)*180/pi  # FOV in deg
# https://www.scantips.com/lights/fieldofview.html
# ISS FOV:  140º diagonal -> 7.87mm FF
# Moon FOV: 1.9º diagonal -> 1300mm FF


# READ WORLD COORDINATES
DT=data.table(map_data("world"))  # long/lat pairs for all countries
DT=DT[, .(num=.N), by=.(long, lat)]  # summarize to deduplicate points

# Deg to rad conversion
DT$phi=DT$long*pi/180
DT$theta=DT$lat*pi/180

# Polar to XYZ coordinates conversion
DT$x=polar2x(Rearth, DT$phi, DT$theta)
DT$y=polar2y(Rearth, DT$phi, DT$theta)
DT$z=polar2z(Rearth, DT$phi, DT$theta)
DT=DT[, list(x, y, z)]  # clean dataframe with 3 columns (x,y,z)


# 1/6: Nations argue...
background=LoadBitmap("background1.png")  # "Nations argue..."
NFRAMES=367
Offset=0
NTURNS=1

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS
    scale=frame/NFRAMES  # global scaling
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot=scaleXYZ(DTplot, sx=scale, sy=scale, sz=scale)
    DTplot$z = DTplot$z + dz  # Earth along Z axis
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points
    
    # Empty bitmap
    img=background*(1-frame/(NFRAMES-1))  # img=NewBitmap(DIMX, DIMY)
    
    # Hidden parts "algorithm":
    #  2. Draw Earth (solid globe)
    #  3. Draw visible maps
    
    # 2. Draw Earth (solid globe)
    img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS*scale,
                       valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)
    
    # 3. Draw visible maps
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points
    
    print(paste0("Part 5/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
    ))
}



# 2/6: Blaming each other...
background=LoadBitmap("background2.png")  # "Blaming each other..."
NFRAMES=350
Offset=367
NTURNS=1

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot$z = DTplot$z + dz  # Earth along Z axis
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points
    
    # Empty bitmap
    img=background*(1-frame/(NFRAMES-1))  # img=NewBitmap(DIMX, DIMY)
    
    # Hidden parts "algorithm":
    #  2. Draw Earth (solid globe)
    #  3. Draw visible maps
    
    # 2. Draw Earth (solid globe)
    img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS,
                       valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)
    
    # 3. Draw visible maps
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points
    
    print(paste0("Part 5/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
    ))
}



# 3/6: And inevitably...
background=LoadBitmap("background3.png")  # "And inevitably..."
NFRAMES=284
Offset=717
NTURNS=1

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot$z = DTplot$z + dz  # Earth along Z axis
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points
    
    # Empty bitmap
    img=background*(1-frame/(NFRAMES-1))  # img=NewBitmap(DIMX, DIMY)
    
    # Hidden parts "algorithm":
    #  2. Draw Earth (solid globe)
    #  3. Draw visible maps

    # 2. Draw Earth (solid globe)
    img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS,
                       valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)
    
    # 3. Draw visible maps
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points

    print(paste0("Part 5/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
    ))
}



# 4/6: STOP
NFRAMES=85
Offset=1001
NTURNS=0

# Rotation and re allocation of maps
DTplot=DT  # initial position of map points
DTplot$z = DTplot$z + dz  # Earth along Z axis

# Distance from each map point to observation point (0,0,0)
DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points

# Empty bitmap
img=NewBitmap(DIMX, DIMY)

# 2. Draw Earth (solid globe)
img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS,
                   valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)

# 3. Draw visible maps
DTplot$factor=f/DTplot$z
DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS

    print(paste0("Part 4/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
    ))
}



# 5/6: War breaks out...
TILT=pi/2  # X tilt for better rendering of North hemisphere
background=LoadBitmap("background5.png")  # "War breaks out..."
NFRAMES=1914
Offset=1086
NTURNS=6

# READ ICBM DATA AND PRECALCULATE ALL TRAJECTORIES
icbm=data.table(read.csv2("icbm.csv"))

# l=Launch, t=Target
icbm$rl=Rearth
icbm$phil=icbm$long_launch*pi/180  # longitude (E/W) in rad
icbm$thetal=icbm$lat_launch*pi/180  # latitude (N/S) in rad

icbm$rt=Rearth
icbm$phit=icbm$long_target*pi/180  # longitude (E/W) in rad
icbm$thetat=icbm$lat_target*pi/180  # latitude (N/S) in rad

icbm$xl=polar2x(icbm$rl, icbm$phil, icbm$thetal)
icbm$yl=polar2y(icbm$rl, icbm$phil, icbm$thetal)
icbm$zl=polar2z(icbm$rl, icbm$phil, icbm$thetal)

icbm$xt=polar2x(icbm$rt, icbm$phit, icbm$thetat)
icbm$yt=polar2y(icbm$rt, icbm$phit, icbm$thetat)
icbm$zt=polar2z(icbm$rt, icbm$phit, icbm$thetat)

# Calculation of arc between Launch and Target using dot product
# (radius of l,t ground locations = Rearth)
icbm$dalpha=acos((icbm$xl*icbm$xt+icbm$yl*icbm$yt+icbm$zl*icbm$zt) / Rearth^2)
maxarc=max(icbm$dalpha)  # longest range missile
icbm$Ni=round(icbm$dalpha/maxarc*NFRAMES)
icbm$Hi=icbm$dalpha/maxarc*MAXH  # same launch angle condition for all ICBM's

# Create list of dataframes (each trajectory = 1 dataframe)
NTRAJ=nrow(icbm)  # number of trajectories
traj=list()  # list of trajectories
Npoints=c()  # vector with total length of each trajectory
for (i in 1:NTRAJ) {
    Npoints[i]=icbm$Ni[i]
    
    # Missile trajectory linear interpolation in XYZ coordinates
    traj[[i]]=data.table(
        'x'=seq(from=icbm$xl[i], to=icbm$xt[i], length.out=Npoints[i]),
        'y'=seq(from=icbm$yl[i], to=icbm$yt[i], length.out=Npoints[i]),
        'z'=seq(from=icbm$zl[i], to=icbm$zt[i], length.out=Npoints[i]))
    
    # Conversion to Polar -> add altitude to r
    traj[[i]]$r=Rearth+icbm$Hi[i]*sin(seq(from=0, to=pi, length.out=Npoints[i]))
    traj[[i]]$phi=xyz2phi(traj[[i]]$x, traj[[i]]$y, traj[[i]]$z)
    traj[[i]]$theta=xyz2theta(traj[[i]]$x, traj[[i]]$y, traj[[i]]$z)
    
    # Final conversion back to XYZ
    traj[[i]]$xp=polar2x(traj[[i]]$r, traj[[i]]$phi, traj[[i]]$theta)
    traj[[i]]$yp=polar2y(traj[[i]]$r, traj[[i]]$phi, traj[[i]]$theta)
    traj[[i]]$zp=polar2z(traj[[i]]$r, traj[[i]]$phi, traj[[i]]$theta)
    
    # Cleaning
    traj[[i]]=traj[[i]][, list(xp, yp, zp)]
    colnames(traj[[i]])=c('x','y','z')
}

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot=rotateX(DTplot, theta=-TILT*frame/(NFRAMES-1))
    DTplot$z = DTplot$z + dz  # Earth along Z axis
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points
    
    trajplot=list()  # keep only points to be plotted in the frame
    for (i in 1:NTRAJ) {
        # Rotation and re allocation of trajectories
        lastpoint=min(frame+1, Npoints[i])  # lastpoint points to plot
        trajplot[[i]]=traj[[i]][1:lastpoint,] # initial position of trajectory points
        trajplot[[i]]=rotateY(trajplot[[i]], theta=-theta)
        trajplot[[i]]=rotateX(trajplot[[i]], theta=-TILT*frame/(NFRAMES-1))
        trajplot[[i]]$z = trajplot[[i]]$z + dz  # Earth along Z axis
        
        # Distance from each trajectory point to observation point (0,0,0)
        trajplot[[i]]$dist=(trajplot[[i]]$x^2+trajplot[[i]]$y^2+trajplot[[i]]$z^2)^0.5
        
        # Grayscale trajectories and identify nukes
        trajplot[[i]]$grayscale=row(trajplot[[i]][,1])/lastpoint  # range 0..1
        trajplot[[i]]$boom=0
        if (lastpoint==Npoints[i]) trajplot[[i]]$boom[Npoints[i]]=1  # nuke
    }
    
    # Empty bitmap
    img=background*(1-frame/(NFRAMES-1))  # img=NewBitmap(DIMX, DIMY)
    
    # Hidden parts "algorithm":
    #  1. Draw trajectories more distant than Earth
    #  2. Draw Earth (solid globe)
    #  3. Draw visible maps
    #  4. Draw trajectories closer than Earth
    
    # 1. Draw trajectories more distant than Earth
    for (i in 1:NTRAJ) {
        trajplottmp=trajplot[[i]][trajplot[[i]]$dist>distmax]  # distant points
        if (nrow(trajplottmp)>0) img=PlotTraj()
    }
    
    # 2. Draw Earth (solid globe)
    img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS,
                       valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)

    # 3. Draw visible maps
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points
    
    # 4. Draw trajectories closer than Earth
    for (i in 1:NTRAJ) {
        trajplottmp=trajplot[[i]][trajplot[[i]]$dist<=distmax]  # close points
        if (nrow(trajplottmp)>0) img=PlotTraj()
    }
    
    print(paste0("Part 5/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
                          ))
}



# 6/6: FADE OUT
TILT=pi/2  # X tilt for better rendering of North hemisphere
NFRAMES=339
Offset=3000
NTURNS=1

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES*NTURNS
    scale=1-frame/(NFRAMES-1)  # global scaling
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot=rotateX(DTplot, theta=-TILT)  # keep pi/2 tilted
    DTplot=scaleXYZ(DTplot, sx=scale, sy=scale, sz=scale)
    DTplot$z = DTplot$z + dz  # Earth along Z axis
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[DTplot$dist<=distmax]  # keep only visible points
    
    trajplot=list()  # keep only points to be plotted in the frame
    for (i in 1:NTRAJ) {
        # Rotation and re allocation of trajectories
        # lastpoint=min(frame+1, Npoints[i])  # lastpoint points to plot
        lastpoint=Npoints[i]  # now all nukes finished
        trajplot[[i]]=traj[[i]][1:lastpoint,] # initial position of trajectory points
        trajplot[[i]]=rotateY(trajplot[[i]], theta=-theta)
        trajplot[[i]]=rotateX(trajplot[[i]], theta=-TILT)  # keep pi/2 tilted
        trajplot[[i]]=scaleXYZ(trajplot[[i]], sx=scale, sy=scale, sz=scale)       
        trajplot[[i]]$z = trajplot[[i]]$z + dz  # Earth along Z axis
        
        # Distance from each trajectory point to observation point (0,0,0)
        trajplot[[i]]$dist=(trajplot[[i]]$x^2+trajplot[[i]]$y^2+trajplot[[i]]$z^2)^0.5
        
        # Grayscale trajectories and identify nukes
        trajplot[[i]]$grayscale=row(trajplot[[i]][,1])/lastpoint  # range 0..1
        trajplot[[i]]$boom=0  # don't plot nukes
        # if (lastpoint==Npoints[i]) trajplot[[i]]$boom[Npoints[i]]=1  # nuke
    }
    
    # Empty bitmap
    img=NewBitmap(DIMX, DIMY)
    
    # Hidden parts "algorithm":
    #  1. Draw trajectories more distant than Earth
    #  2. Draw Earth (solid globe)
    #  3. Draw visible maps
    #  4. Draw trajectories closer than Earth
    
    # 1. Draw trajectories more distant than Earth
    for (i in 1:NTRAJ) {
        trajplottmp=trajplot[[i]][trajplot[[i]]$dist>distmax]  # distant points
        if (nrow(trajplottmp)>0) img=PlotTraj()
    }
    
    # 2. Draw Earth (solid globe)
    img=DrawGradCircle(img, NCOLDIV2, NROWDIV2, RADIUS*scale,
                       valmin=GRAYGLOBEMIN, valmax=GRAYGLOBEMAX)
    
    # 3. Draw visible maps
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=GRAYMAP  # draw points
    
    # 4. Draw trajectories closer than Earth
    for (i in 1:NTRAJ) {
        trajplottmp=trajplot[[i]][trajplot[[i]]$dist<=distmax]  # close points
        if (nrow(trajplottmp)>0) img=PlotTraj()
    }
    
    print(paste0("Part 6/6: ", frame, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
    
    SaveBitmap(img, paste0("img", ifelse(frame+Offset<10, "000",
                                  ifelse(frame+Offset<100, "00",
                                  ifelse(frame+Offset<1000, "0", ""))),
                                  frame+Offset
    ))
}


# Building MP4:
# ffmpeg -framerate 30 -i img%4d.png -i AndOnePlayingDead.wav
#        -c:v libx264 -crf 18 -pix_fmt yuv420p oppenheimer.mp4



################################################################################

# ICBM TRAJECTORY APPROXIMATION

DIMX=1920  # Full HD: 1920 x 1080 pixels
DIMY=1080
x0=200
y0=DIMY/2
R=200
img=NewBitmap(DIMX, DIMY)
img=DrawCircle(img, x0, y0, r=R, val=0.5, fill=TRUE)
img=DrawLine(img, x0, 1, x0, DIMY, val=0.5)
img=DrawLine(img, 1, y0, DIMX, y0, val=0.5)

HMISSILE=100
FACTOR=1  # 1.5
N=200
for (dtheta in seq(from=pi/4/3, to=pi, length.out=12)) {
    for (t in 0:N) {
        theta=dtheta*t/N
        h=HMISSILE*dtheta/pi * abs(sin(theta/dtheta*pi))^FACTOR  # abs() avoids
        x=(R+h)*sin(theta)  # rounding errors
        y=(R+h)*cos(theta)
        img=DrawPoint(img, x+x0, y+y0)
    }
}

ShowBitmap(img)
SaveBitmap(img, "icbmtrajectory.png")
