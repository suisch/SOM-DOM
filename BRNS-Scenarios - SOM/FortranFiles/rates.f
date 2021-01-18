c      
c     SUBROUTINE rates
c      
      subroutine rates(j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            r(1,j) = miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)
            r(2,j) = miuMax_F1*sp(1,j)/(k_sF1+sp(1,j))*sp(6,j)
            r(3,j) = miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))*sp(7,j)
            r(4,j) = miuMax_AE*sp(2,j)/(k_sAE+sp(2,j))*sp(5,j)
            r(5,j) = -b_AE*sp(5,j)
            r(6,j) = -b_F1*sp(6,j)
            r(7,j) = -b_F2*sp(7,j)
      end
