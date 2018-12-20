            subroutine potential_h3(r,en)
                    implicit none
                    integer :: i
                    integer, parameter :: nn=3

                    real*8, intent(in) :: r(nn)
                    real*8, intent(out) :: en

                    real*8 :: dv(nn),rr(nn),dva(nn),da(nn),dq(nn),dvt(nn)
                    real*8 :: vt(nn),vs(nn),a(nn),q(nn),dvs(nn)
                    real*8 :: ex7,vv,fsw,vlr

                    real*8 :: beta1,beta2,beta3
                    data beta1,beta2,beta3 /0.52,0.052,0.79/
                    real*8 :: a2,a3,a4,a5
                    data a2,a3,a4,a5 /0.0012646477,-0.0001585792,0.0000079707, &
                                  & -0.0000001151/
                    real*8 :: b11,b12,b22,b23
                    data b11,b12,b22,b23 /3.0231771503,-1.08935219,&
                                  & 1.7732141742,-2.0979468223/
                    real*8 :: b24,b31,b32,b41,b42
                    data b24,b31,b32,b41,b42 /-3.978850217,0.4908116374,&
                                  & -0.8718696387,0.1612118092,-0.1273731045/
                    real*8 :: b51,b52
                    data b51,b52 /-13.3599568553,0.9877930913/
                    real*8 :: alf5
                    data alf5 /0.0035/
                    real*8 :: c6,c8
                    data c6,c8 /6.89992032,219.9997304/
                    real*8 :: xx(8)
                    data xx /2.80465,2.86869,1.93426,-0.85360,0.24768, &
                         & -0.02013,2.80521,0.01927/

                    rr=r-1.40059

                    do i=1,nn

                       ex7=exp(-xx(7)*rr(i))
                       vv=rr(i)*xx(8)+xx(6)
                       vv=rr(i)*vv+xx(5)
                       vv=rr(i)*vv-xx(4)
                       vv=rr(i)*vv+xx(3)
                       vv=rr(i)*vv+xx(2)
                       vv=rr(i)*vv+xx(1)
                       vv=rr(i)*vv+1.0
                       vs(i)=-0.174475*ex7*vv


                       vv=rr(i)*7.0*xx(8)+6.0*xx(6)
                       vv=rr(i)*vv+5.0*xx(5)
                       vv=rr(i)*vv-4.0*xx(4)
                       vv=rr(i)*vv+3.0*xx(3)
                       vv=rr(i)*vv+2.0*xx(2)
                       vv=rr(i)*vv+xx(1)
                       dvs(i)=-0.174475*ex7*vv-xx(7)*vs(i)

                       fsw=exp(-0.011*(r(i)-10.0)**4)
                       if (r(i).gt.10.0) fsw=1.0

                       dfsw=-0.044*(r(i)-10.0D0)**3*fsw
                       if (r(i).gt.10.0) dfsw=0.0

                       vlr=-c6/r(i)**6-c8/r(i)**8



                    enddo



                    





                    






            end subroutine potential_h3

