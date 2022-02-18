rm(list=ls())
library(ggplot2)
library(plotly)
################################################################################
#### Groundwater modeling
#### 2-D Transient Model
#### unconfined aquifer
#### Solved with Successive Over Relaxation scheme
#########################################################################
K = 5e-3                                          # hydraulic conductivity, m/s
dx = 100                                          # grid length, m
dy = 100                                          # grid width, m
nx = 21                                           # grids in x direction
ny = 16                                           # grids in y direction
theta = 0.4                                       # porosity

zbot = -100                                       # grid bottom elevation (aquifer depth), m
ztop = 0                                          # grid top elevation (aquifer height), m
theta = 0.4                                       # porosity

tt <- 365*24*3600*500                             # simulation time, 500 years in seconds
dt <- 3600             

nt <- ceiling(tt/dt)
dt <- tt/nt
t <- seq(0,tt,dt)
Hl= 95                                            # water head at west side
Hr= 80                                            # river water surface elevation

Sy = 0.25                                         # specific yield
S0 = Sy                                           # because of unconfined aquifer, S0 = Sy

#### source term well pumping rate
W = matrix(0,nx,ny)                               # source term
W[9:13,8] <- -0.91                                # only 2 dimensions because of a fixed pumping rate
# pumping location and pumping rate is 0.05
Ax = K*dy/dx                                      # conductance across the entire grid
Ay = K*dx/dy

Aw <- matrix(Ax,nx,ny)                            # conductance at each cell face
Ae <- matrix(Ax,nx,ny)
As <- matrix(Ay,nx,ny)
An <- matrix(Ay,nx,ny)
Aw[,1] <- 2*Ax                                    # mult by 2, right at cell face
Ae[,ny] <- 2*Ax
As[1,] <- 2*Ay
An[nx,] <- 2*Ay


H = matrix(95,nx,ny)                 # initial water head
Hold=H                               # store water head at previous iteration
qx = matrix(0,nx,ny+1)               # flow at each grid interface in x direction
qy = matrix(0,nx+1,ny)               # flow at each grid interface in y direction


#############################################
#### define contaminant properties and initial condition
#### use mercury as an example
C0 = 1200                            # initial porewater concentration, ug/L
C_RG = 3.8e-2                        # Fresh water quality standard, ug/L
Disp = 6.54e-10                      # molecular diffusion, m2/s
rho = 1.56                           # soil bulk density, kg/m3
Kd = 10^2.2                          # PAH-soil linear partition coefficient, L/kg
R = 1+rho*Kd/theta                   # retardation factor
Cpw = matrix(0,nx,ny)


epsilon = 1e-3
max.iter = 100
omega = 1.5                          # coefficient controls the over-relaxation
H_out <- list()
qx_out <- list()
qy_out <- list()
Cpw_out <- list()

t.save = 0
ind.save = 0
#### assign a new matrix with dimension (nx+2) * (ny+2)
#### boundary conditions are given in the edge cells
Z <- matrix(95,nx+2,ny+2)
Zt <- Z
CZ <- matrix(0,nx+2,ny+2)


#### start simulation
for (k in 1:nt) {
  t.save = t.save+1
  Z[2:(nx+1),2:(ny+1)] = H                        # assign water head to internal grids
  Z[,1]=Hl
  Z[,ny+2]=Hr
  Z[1,2:(ny+1)]=H[1,]
  Z[nx+2,2:(ny+1)]=H[nx,]
  Zt = Z                         # store the water head at previous time step
  #### fixed initial concentration for grids represent dump site
  Cpw[10:12,6] <- C0                               # this represents the cells of the disposal site
  CZ[2:(nx+1),2:(ny+1)] = Cpw
  
  err = 100                                        
  niter = 0                                        # initialize iteration
  beta = S0*dx*dy/dt
  while (err > epsilon & niter <= max.iter) {
    Z[2:(nx+1),2:(ny+1)] = H
    Hold = H                                           # water head at current iteration
    for (i in 2:(nx+1)) {
      for (j in 2:(ny+1)) {
        Cw = Aw[i-1,j-1]*(Z[i,j-1]+Z[i,j])/2
        Ce = Ae[i-1,j-1]*(Z[i,j]+Z[i,j+1])/2
        Cs = As[i-1,j-1]*(Z[i-1,j]+Z[i,j])/2
        Cn = An[i-1,j-1]*(Z[i,j]+Z[i+1,j])/2
        Cp = beta + Cw + Ce + Cs + Cn                 # part of LHS 
        Z[i,j] <- (1-omega)*Z[i,j]+omega/Cp*(beta*Zt[i,j]+Cw*Z[i,j-1]+Ce*Z[i,j+1]+Cs*Z[i-1,j]+Cn*Z[i+1,j]+W[i-1,j-1])
      }
    }
    H = Z[2:(nx+1),2:(ny+1)]                          # get new water head
    err = sum(sum(abs(H-Hold)))                       # update error
    niter = niter + 1                                 # iteration is increased by one
  }
  if (t.save==24) {
    t.save=0
    ind.save=ind.save+1
    H_out[[ind.save]] = H
    
    #### calculate flow
    qx[,1] <- Aw[,1]*(Hl+H[,1])/2*(Hl-H[,1])
    Hbar <- (H[,1:(ny-1)]+H[,2:ny])/2
    qx[,2:ny] <- Aw[,2:ny]*Hbar*(H[,1:(ny-1)]-H[,2:ny])
    qx[,ny+1] <- Ae[,ny]*(H[,ny]+Hr)/2*(H[,ny]-Hr)
    
    qy[1,] <- 0
    Hbar <- (H[1:(nx-1),]+H[2:nx,])/2
    qy[2:nx,] <- As[2:nx,]*Hbar*(H[1:(nx-1),]-H[2:nx,])
    qy[nx+1,] <- 0
    qx_out[[ind.save]] <- qx
    qy_out[[ind.save]] <- qy
  }
  ##### start the contaminant transport simulation
  ##### use explicit scheme
  ##### check time step
  # dtd <- dx*dx/Disp/4247222222.22
  # vx <- qx[1:nx,1:ny]/H/dy
  # vy <- qy[1:nx,1:ny]/H/dx
  # dtvx <- dx/max(vx)
  # dtvy <- dy/max(vy)
  
  dv = dx * dy * Z
  Rdt = (R * dv) / dt
  Wc = W[i-1,j-1] * CZ[i,j]                     
  dispx = (dy * Disp * Z) / dx                # dispersion term in x dimension
  dispy = (dx * Disp * Z) / dy 
  
  for(i in 2:(nx+1)) {
    for(j in 2:(ny+1)) {
      #CZ[i,j] = C0                             # C must be initialized?
      CZ[i,j] = abs(((Rdt[i,j]*CZ[i,j])-(CZ[i,j]*(qx[i-1,j-1]+qy[i-1,j-1]+(2*dispx[i,j])+(2*dispy[i,j])))
                 +(qx[i-1,j-1]*(CZ[i,j-1]))
                 +(qy[i-1,j-1]*(CZ[i-1,j]))
                 +(dispx[i,j]*(CZ[i,j-1]+CZ[i,j+1]))
                 +(dispy[i,j]*(CZ[i-1,j]+CZ[i+1,j]))
                 +Wc)/Rdt[i,j])
    }
  }
  Cpw <- CZ[2:(nx+1), 2:(ny+1)]
  
  if (t.save==365) {                                             # run this after the first for loop
    t.save=0
    ind.save=ind.save+1
    H_out[[ind.save]] = H
    Cpw_out[[ind.save]] = Cpw
    qx_out[[ind.save]] <- qx
    qy_out[[ind.save]] <- qy
  }
}  


X <- seq(0.5*dx,(nx-0.5)*dx,dx)       # length along x direction
Y <- seq(0.5*dy,(ny-0.5)*dy,dy)
plot_ly(x = X, y = Y, z = H, type = "contour", contours = list(showlabels = TRUE))
plot_ly(x = X, y = Y, z = Cpw, type = "contour", contours = list(showlabels = TRUE))

