        program sho_dynamics
                implicit none
                integer :: i,nst
                real*8 :: m,k
                real*8 :: t_init,t_final,x0,p0
                real*8 :: x,p,t,dt 
                real*8 :: xk1,xk2,xk3,xk4
                real*8 :: pk1,pk2,pk3,pk4
                common/const/m,k

                open(1,file="input.dat",status="old")
                open(10,file="output.dat",status="unknown")

                read(1,*)m
                read(1,*)k
                read(1,*)dt
                read(1,*)x0,p0

                t_init=0.0
                t_final=4.0*(4.0*atan(1.0)*sqrt(m/k))
                p0=m*p0

                nst=int((t_final-t_init)/dt)

                x=x0
                p=p0
                t=t_init

                write(10,1)t,x,p

                do i=1,nst

                   call fex(x,p,xk1,pk1)
                   call fex(x+dt*xk1/2.0,p+dt*pk1/2.0,xk2,pk2)
                   call fex(x+dt*xk2/2.0,p+dt*pk2/2.0,xk3,pk3)
                   call fex(x+dt*xk3,p+dt*pk3,xk4,pk4)

                   x=x+dt*(xk1+2.0*xk2+2.0*xk3+xk4)/6.0
                   p=p+dt*(pk1+2.0*pk2+2.0*pk3+pk4)/6.0

                   t=t+dt

                   write(10,1)t,x,p

                enddo

                close(1)
                close(2)

 1              format(3f12.6)
                stop

        end program sho_dynamics

        subroutine fex(x,p,xdot,pdot)        

                implicit none
                real*8 :: m,k
                real*8, intent(in) :: x,p
                real*8, intent(out) :: xdot,pdot
                common/const/m,k

                xdot=p/m
                pdot=-k*x

                return
        end subroutine
