MODULE smat

USE tdsim

REAL, ALLOCATABLE :: S_matrix(:,:,:),xx(:,:,:),yy(:,:,:),zz(:,:,:),F_matrix(:,:,:),c(:)
REAL :: T,INTG=0.0

CONTAINS 
SUBROUTINE matgen(c,xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,numvar)	
REAL, ALLOCATABLE :: c(:)

allocate(F_matrix(numvar,numvar,numvar))
IF(allocated(F_matrix) .NEQV. .TRUE.) THEN
	write(*,*) "Unable to allocate memory"
END IF
allocate(S_matrix(numvar,numvar,numvar))
IF(allocated(S_matrix) .NEQV. .TRUE.) THEN
	write(*,*) "Unable to allocate memory"
END IF
allocate(xx(numvar,numvar,numvar))
IF(allocated(xx) .NEQV. .TRUE.) THEN
	write(*,*) "Unable to allocate memory"
END IF
allocate(yy(numvar,numvar,numvar))
IF(allocated(yy) .NEQV. .TRUE.) THEN
	write(*,*) "Unable to allocate memory"
END IF
allocate(zz(numvar,numvar,numvar))
IF(allocated(zz) .NEQV. .TRUE.) THEN
	write(*,*) "Unable to allocate memory"
END IF

	DO l=1,numvar
		DO j=1,numvar
			DO k=1,numvar
				S_matrix(l,j,k)=c(l)*c(j)*c(k)
			END DO
		END DO
	END DO
  
	DO l=1,(numvar)
		DO j=1,(numvar)
			DO k=1,(numvar)
				xx(l,j,k)= xmin+(k-1)*dx
			END DO
		END DO
	END DO

	DO l=1,(numvar)
		DO k=1,(numvar)
			DO j=1,(numvar)
				yy(l,j,k)= ymin+(j-1)*dy
			END DO
		END DO
	END DO

	DO j=1,(numvar)
		DO k=1,(numvar)
			DO l=1,(numvar)
				zz(l,j,k)= zmin+(l-1)*dz
			END DO
		END DO
	END DO

	DO l=1,(numvar)
		DO j=1,(numvar)
			DO k=1,(numvar)
				CALL fgen(xx(l,j,k),yy(l,j,k),zz(l,j,k),T)
				F_matrix(l,j,k) = T
				INTG = INTG + S_matrix(l,j,k)*F_matrix(l,j,k)
			END DO 
		END DO
	END DO
CALL integ(INTG)
END SUBROUTINE matgen

SUBROUTINE fgen(x,y,z,r)
!Define function to be utilised
REAL :: x,y,z,r
r = x**3 - 3*y*z
END SUBROUTINE fgen

SUBROUTINE integ(INTG)
REAL :: INTG
INTG=((dx*dy*dz)/27)*(INTG)
RETURN
END SUBROUTINE integ

END MODULE smat


