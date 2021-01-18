c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = -1/Y_AE*miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)-miuM
     +ax_F1*sp(1,j)/(k_sF1+sp(1,j))*sp(6,j)/Y_F1-miuMax_F2*sp(1,j)/(k_sF
     +2+sp(1,j))*sp(7,j)/Y_F2
            drdc = -1/Y_AE*miuMax_AE/(k_sAE+sp(1,j))*sp(5,j)+1/Y_AE*miuM
     +ax_AE*sp(1,j)/(k_sAE+sp(1,j))**2*sp(5,j)-miuMax_F1/(k_sF1+sp(1,j))
     +*sp(6,j)/Y_F1+miuMax_F1*sp(1,j)/(k_sF1+sp(1,j))**2*sp(6,j)/Y_F1-mi
     +uMax_F2/(k_sF2+sp(1,j))*sp(7,j)/Y_F2+miuMax_F2*sp(1,j)/(k_sF2+sp(1
     +,j))**2*sp(7,j)/Y_F2
            endif
            if (isp.eq.2) then
            rat = -1/Y_AE*miuMax_AE*sp(2,j)/(k_sAE+sp(2,j))*sp(5,j)
            drdc = -1/Y_AE*miuMax_AE/(k_sAE+sp(2,j))*sp(5,j)+1/Y_AE*miuM
     +ax_AE*sp(2,j)/(k_sAE+sp(2,j))**2*sp(5,j)
            endif
            if (isp.eq.3) then
            rat = -2*miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)*(-Y_AE*al
     +pha-2*Y_AE+1)/Y_AE
            drdc = 0
            endif
            if (isp.eq.4) then
            rat = (1-Y_AE)*miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)/Y_A
     +E+(1-Y_F1-(-Y_F1*gammaB+gammaGlu)/gammaEth)*miuMax_F1*sp(1,j)/(k_s
     +F1+sp(1,j))*sp(6,j)/Y_F1+(-Y_F2*gammaB+gammaGlu)/gammaLac*miuMax_F
     +2*sp(1,j)/(k_sF2+sp(1,j))*sp(7,j)/Y_F2
            drdc = 0
            endif
            if (isp.eq.5) then
            rat = miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)-b_AE*sp(5,j)
            drdc = miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))-b_AE
            endif
            if (isp.eq.6) then
            rat = miuMax_F1*sp(1,j)/(k_sF1+sp(1,j))*sp(6,j)-b_F1*sp(6,j)
            drdc = miuMax_F1*sp(1,j)/(k_sF1+sp(1,j))-b_F1
            endif
            if (isp.eq.7) then
            rat = miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))*sp(7,j)-b_F2*sp(7,j)
            drdc = miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))-b_F2
            endif
      end
