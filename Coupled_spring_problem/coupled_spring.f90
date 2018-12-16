        program Coupled_spring_problem
                implicit none
                integer, parameter :: n=2
                integer :: i,j,nst
                real*8 :: m(n),k(n),del(n),mu(n)
                real*8 :: x0(n),v0(n)
                real*8 :: x(n),v(n)
                real*8 :: t_init,t_final,t,dt
                real*8 :: xk1(n),xk2(n),xk3(n),xk4(n)
                real*8 :: vk1(n),vk2(n),vk3(n),vk4(n)
                common/const/m,k,del,mu

                open(1,file="input.dat",status="old")
                open(10,file="output.dat",status="unknown")

                t_init=0.0
                
                read(1,*)(m(i),i=1,n)
                read(1,*)(k(i),i=1,n)
                read(1,*)(x0(i),i=1,n)
                read(1,*)(v0(i),i=1,n)
                read(1,*)(del(i),i=1,n)
                read(1,*)(mu(i),i=1,n)
                read(1,*)t_final,dt

                nst=int((t_final-t_init)/dt)

                x=x0
                v=v0
                t=t_init

                write(10,1)t,(x(i),v(i),i=1,n)

                do i=1,nst

                   call fex(x,v,xk1,vk1)
                   call fex(x+dt*xk1/2.0,v+dt*vk1/2.0,xk2,vk2)
                   call fex(x+dt*xk2/2.0,v+dt*vk2/2.0,xk3,vk3)
                   call fex(x+dt*xk3,v+dt*vk3,xk4,vk4)

                   x=x+dt*(xk1+2.0*xk2+2.0*xk3+xk4)/6.0
                   v=v+dt*(vk1+2.0*vk2+2.0*vk3+vk4)/6.0

                   t=t+dt

                   write(10,1)t,(x(j),v(j),j=1,n)

                enddo

                close(1)
                close(2)

 1              format(5f12.6)
                stop

        end program Coupled_spring_problem

        subroutine fex(x,v,xdot,vdot)        

                implicit none
                integer, parameter :: n=2
                real*8 :: m(n),k(n),del(n),mu(n)
                real*8, intent(in) :: x(n),v(n)
                real*8, intent(out) :: xdot(n),vdot(n)
                common/const/m,k,del,mu

                xdot(1)=v(1)
                vdot(1)=-k(1)*x(1)/m(1)-k(2)*(x(1)-x(2))/m(1)
                vdot(1)=vdot(1)-del(1)*v(1)
                vdot(1)=vdot(1)+mu(1)*x(1)**3+mu(2)*(x(1)-x(2))**3
                xdot(2)=v(2)
                vdot(2)=-k(2)*(x(2)-x(1))/m(2)
                vdot(2)=vdot(2)-del(2)*v(2)
                vdot(2)=vdot(2)+mu(2)*(x(2)-x(1))**3

                return
        end subroutine
