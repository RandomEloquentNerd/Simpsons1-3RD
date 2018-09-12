!Program to calculate the integral of a given function using the Simpson's 1/3rd rule
!I=(h/3)(f(x0)+f(xn)+4*sum(f(xodd))+2*sum(f(xeven)))
!	author: Kushan
!	date: 14/08/18
PROGRAM simpsons

USE tdsim
USE smat

REAL, allocatable :: func(:)
INTEGER :: m,stat1,stat2,k=1,p
REAL :: div=0.0,I=0.0
CHARACTER :: user,data_set

write(*,*) "Multi-Integral Simpsons or just the usual, my liege? (m,u)"
read(*,*) user

!DO i=1,m
!write(*,*) "Enter the value for the variable x(" ,i, ")now:"
!read(*,*) x(i)
!write(*,*) "Enter the value for the function at x(" ,i, ")now:"
!read(*,*) func(i)
!END DO
SELECT CASE(user)
	CASE('u')
		write(*,*) "Which data set do you wish to tap into (1/2)?"
		read(*,*) data_set
		SELECT CASE(data_set)
			CASE('1')
				write(*,*) "Enter the number of values in your discrete distribution:"
				read(*,*) m
				allocate(x(m))
				allocate(func(m))
				IF(allocated(x) .NEQV. .TRUE.) THEN
					write(*,*) "Unable to allocate memory"
				END IF
				IF(allocated(func) .NEQV. .TRUE.) THEN
					write(*,*) "Unable to allocate memory"
				END IF
				DO
					IF(stat1/=0) EXIT
					open (unit=13,file='ipxa',ACTION='READ')
					read(13,*,IOSTAT=stat1) x(k)	
					open (unit=97,file='ipf(x)a',ACTION='READ')	
					read(97,*,IOSTAT=stat1) func(k)
					k=k+1
				END DO

				div=(x(m)-x(1))/(m-1)
				I=(div/3.0)*(func(1)+func(m))

				DO j=2,m-1
					IF(modulo(j,2) .EQ. 0.0) THEN
						I=I+4.0*div*func(j)/3.0
					ELSE
						I=I+2.0*div*func(j)/3.0
					END IF
				END DO
				open(unit=96,file='opI',POSITION='APPEND',ACTION='READWRITE')
				write(96,*) '1',I
				deallocate(x)
				deallocate(func)
				close(13)
				close(97)
				close(96)
			CASE('2')
				write(*,*) "Enter the number of values in your discrete distribution:"
				read(*,*) m
				allocate(x(m))
				allocate(func(m))
				IF(allocated(x) .NEQV. .TRUE.) THEN
					write(*,*) "Unable to allocate memory"
				END IF
				IF(allocated(func) .NEQV. .TRUE.) THEN
					write(*,*) "Unable to allocate memory"
				END IF
				k=1
				DO
					IF(stat2/=0) EXIT
					open (unit=14,file='ipxb',ACTION='READ')
					read(14,*,IOSTAT=stat2) x(k)	
					open (unit=98,file='ipf(x)b',ACTION='READ')	
					read(98,*,IOSTAT=stat2) func(k)
					k=k+1
				END DO
				
				div=(x(m)-x(1))/(m-1)
				
				I=(div/3.0)*(func(1)+func(m))
				
				DO j=2,m-1
					IF(modulo(j,2) .EQ. 0.0) THEN
						I=I+4.0*div*func(j)/3.0
					ELSE
						I=I+2.0*div*func(j)/3.0
					END IF
				END DO
				
				open(unit=96,file='opI',POSITION='APPEND',ACTION='READWRITE')
				write(96,*) '2',I
				deallocate(x)
				deallocate(func)
				close(14)
				close(98)
				close(96)
			CASE DEFAULT
				write(*,*) "Erroneous input re-run the code"
		END SELECT
	CASE('m')

		write(*,*) "Enter the number of variables you're dealing with:(<3)"
		read(*,*) varnum
		allocate(varmin(varnum))
		IF(allocated(varmin) .NEQV. .TRUE.) THEN
				write(*,*) "Unable to allocate memory"
		END IF
		allocate(varmax(varnum))
		IF(allocated(varmax) .NEQV. .TRUE.) THEN
				write(*,*) "Unable to allocate memory"
		END IF
		write(*,*) "Enter the # of partitions for the variables:"
		read(*,*) partnum
		
		k=1		

		DO
			IF(k-1==varnum) EXIT
			IF(stat/=0) EXIT
				open (unit=10,file='ip',ACTION='READ')
				read(10,*,IOSTAT=stat) varmin(k),varmax(k) 
				write(*,*) varmin(k),varmax(k)				
				k=k+1
		END DO

		xmax=varmax(1);xmin=varmin(1);ymax=varmax(2);ymin=varmin(2);zmax=varmax(3);zmin=varmin(3)

		CALL step(xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,partnum)
		write(*,*) dx,dy,dz
		numvar=partnum+1
		allocate(c(numvar))
		IF(allocated(c) .NEQV. .TRUE.) THEN
			write(*,*) "Unable to allocate memory"
		END IF
		DO k=1,numvar
			IF((k .EQ. 1) .OR. (k .EQ. (numvar))) THEN
				c(k)=1 
			ELSE IF(modulo(k,2) .EQ. 0) THEN
				c(k)=4
			ELSE 
				c(k)=2
			END IF
		END DO
		
		CALL matgen(c,xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,numvar)	

		open(unit=96,file='opI',POSITION='APPEND',ACTION='READWRITE')
		write(96,*) '3',INTG	

	CASE DEFAULT
		write(*,*) "Erroneous value entered, please re-run the code and re-enter one of the designated values for Integrals"
END SELECT
END PROGRAM simpsons
