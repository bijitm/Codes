            subroutine potential_h3(r,en)
                    implicit none
                    integer, parameter :: nn=3
                    real*8, intent(in) :: r(nn)
                    real*8, intent(out) :: en
                    real*8 :: dv(nn),rr(nn),dva(nn),da(nn),dq(nn),dvt(nn)
                    real*8 :: vt(nn),vs(nn),a(nn),q(nn),dvs(nn),beta(nn)
                    data beta(:) /0.52D0,0.052D0,0.79D0/

                    






            end subroutine potential_h3

