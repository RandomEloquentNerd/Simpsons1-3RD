MODULE tdsim

REAL, allocatable :: x(:),y(:),z(:),varmax(:),varmin(:)
REAL :: xmax,xmin,ymax,ymin,zmax,zmin,dx=0.0,dy=0.0,dz=0.0,partnum=0.0
INTEGER :: k1=1,k2=1,k3=1,stat,varnum,numvar=0

CONTAINS 
SUBROUTINE step(xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,partnum)
dx=(xmax-xmin)/(partnum)
dy=(ymax-ymin)/(partnum)
dz=(zmax-zmin)/(partnum)
RETURN
END SUBROUTINE step

END MODULE tdsim


