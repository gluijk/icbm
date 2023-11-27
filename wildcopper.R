# Wild Copper demo tribute
# www.overfitting.net
# https://www.overfitting.net/2023/08/oppenheimer-i-am-become-death-destroyer.html


library(ggmap)  # map_data()
library(data.table)
library(png)


# BITMAP DRAWING FUNCTIONS

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}


# COORDINATES CONVERSION

# Polar to XYZ
polar2x = function(r, phi, theta) { r*cos(theta)*sin(phi)}
polar2y = function(r, phi, theta) { r*sin(theta)}
polar2z = function(r, phi, theta) {-r*cos(theta)*cos(phi)}

rotateY = function(df, theta=0) {  # rotation around Y axis
    df$xr = cos(theta)*df[,1] + sin(theta)*df[,3]
    df$yr = df[,2]    
    df$zr =-sin(theta)*df[,1] + cos(theta)*df[,3]
    df=df[,4:6]  # keep only rotated coords
    colnames(df)=c('x','y','z')
    return(df)
}


################################################################################


# DRAW EARTH

# PHYSICAL PARAMETERS
Rearth=6371.23  # Earth average radius (km)
dz=Rearth*4  # observation point to Earth centre distance (km)
thetamax=acos(Rearth/dz)


# ANIMATION PARAMETERS
DIMX=1920  # Full HD animation
DIMY=1080  # 1920 x 1080 pixels
NCOLDIV2=round(DIMX/2)
NROWDIV2=round(DIMY/2)
TH=0.65  # border between Earth globe and image limits
COLOURGAP=1.2

# Calculate focal length to fit Earth in the final image and FOV
f=min(NCOLDIV2,NROWDIV2)*TH*
    (dz-Rearth*cos(thetamax))/(Rearth*sin(thetamax))
FOV=(pi-2*thetamax)*180/pi  # FOV in deg


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


# ANIMATION
NFRAMES=235  # 54,896s audio track at 30fps
NTURNS=7
background=readPNG("background.png")  # blue plane (Photoshop)

for (frame in 0:(NFRAMES-1)) {
    theta=2*pi*frame/NFRAMES  #*NTURNS
    
    # Rotation and re allocation of maps
    DTplot=DT  # initial position of map points
    DTplot=rotateY(DTplot, theta=-theta)
    DTplot$z = DTplot$z + dz  # Earth Z axis (distance) shift
    DTplot$y = DTplot$y + Rearth/3  # Earth Y axis (vertical) shift
    
    # Distance from each map point to observation point (0,0,0)
    DTplot$dist=(DTplot$x^2+DTplot$y^2+DTplot$z^2)^0.5
    DTplot=DTplot[order(-dist)]  # order points by distance to observer
    MAX=max(DTplot$dist)
    MIN=min(DTplot$dist)

    # Empty bitmap
    img=NewBitmap(DIMX, DIMY)

    # Draw globe
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2  # 3D to 2D projection
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    for (i in 1:nrow(DTplot)) img[round(DTplot$xp[i]),
        round(DTplot$yp[i])]=(1-(DTplot$dist[i]-MIN)/(MAX-MIN))/COLOURGAP+(1-1/COLOURGAP)

    # Draw shadow
    DTplot$yp=-Rearth*1.1*DTplot$factor + NROWDIV2  # projection plane for shadows
    for (i in 1:nrow(DTplot)) img[round(DTplot$xp[i]), round(DTplot$yp[i])]=-1
    
    # Generate colour image
    imgout=background
    img=t(img[,ncol(img):1])
    
    indices=which(img>0)  # globe
    imgout[,,1][indices]=img[indices]
    imgout[,,2][indices]=0
    imgout[,,3][indices]=0
    
    indices=which(img==-1)  # shadows
    imgout[,,1][indices]=0
    imgout[,,2][indices]=0
    imgout[,,3][indices]=0  
    
    for (i in 0:(NTURNS-1)) {
    writePNG(imgout, paste0("img", ifelse(frame+i*NFRAMES<10, "000",
        ifelse(frame+i*NFRAMES<100, "00",
        ifelse(frame+i*NFRAMES<1000, "0", ""))),
        frame+i*NFRAMES, ".png"))
    }

    print(paste0(frame+1, "/", NFRAMES,
                 ", theta=", round(theta*180/pi), "º, ",
                 nrow(DTplot), " points"))
}



# Building MP4:
# ffmpeg -framerate 30 -i img%4d.png -i wildcopper.wav
#        -c:v libx264 -crf 18 -pix_fmt yuv420p wildcopper.mp4
