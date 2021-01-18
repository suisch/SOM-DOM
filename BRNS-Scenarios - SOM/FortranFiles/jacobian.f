c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(1,1) = -1/Y_AE*miuMax_AE/(k_sAE+sp(1,j))*sp(5,j)+1/Y_AE*miuM
     +ax_AE*sp(1,j)/(k_sAE+sp(1,j))**2*sp(5,j)
         pd(2,6) = 0
         pd(1,7) = 0
         pd(2,5) = -1/Y_AE*miuMax_AE*sp(2,j)/(k_sAE+sp(2,j))
         pd(4,5) = 0
         pd(7,1) = -1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*
     +gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gamma
     +Lac)*Y_F2*gammaLac*(Y_F1*gammaB-Y_F1*gammaEth+gammaEth-gammaGlu)/d
     +elt
         pd(7,5) = 0
         pd(3,3) = (Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaLac)*gammaEt
     +h/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEt
     +h-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)/delt
         pd(4,4) = -1/delt
         pd(3,4) = 2*(Y_AE*alpha+2*Y_AE-1)/(Y_F1*gammaB*gammaLac-Y_F1*ga
     +mmaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*ga
     +mmaLac-gammaGlu*gammaLac)*gammaEth*gammaLac/delt
         pd(3,1) = -2*(Y_AE*alpha+2*Y_AE-1)/Y_F1*miuMax_F1/(k_sF1+sp(1,j
     +))*sp(6,j)+2*(Y_AE*alpha+2*Y_AE-1)/Y_F1*miuMax_F1*sp(1,j)/(k_sF1+s
     +p(1,j))**2*sp(6,j)-2*gammaEth*(Y_AE*Y_F2*alpha*gammaB+2*Y_AE*Y_F2*
     +gammaB-Y_AE*alpha*gammaGlu-2*Y_AE*gammaGlu-Y_F2*gammaB+gammaGlu)/(
     +Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-g
     +ammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)/delt
         pd(6,7) = 0
         pd(7,2) = 0
         pd(3,5) = 0
         pd(4,1) = -(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*ga
     +mmaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLa
     +c)/gammaEth/Y_F2/gammaLac*miuMax_F2/(k_sF2+sp(1,j))*sp(7,j)+(Y_F1*
     +gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaE
     +th*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)/gammaEth/Y_F2/gam
     +maLac*miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))**2*sp(7,j)-(Y_F1*gammaB-Y_
     +F1*gammaEth+gammaEth-gammaGlu)/gammaEth/delt
         pd(7,4) = -1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*
     +gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gamma
     +Lac)*gammaEth*Y_F2*gammaLac/delt
         pd(1,4) = 0
         pd(4,7) = -(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*ga
     +mmaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLa
     +c)/gammaEth/Y_F2/gammaLac*miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))
         pd(5,5) = -b_AE-1/delt
         pd(6,6) = -b_F1-1/delt
         pd(7,3) = -1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*
     +gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gamma
     +Lac)*Y_F2*gammaLac*(Y_AE*gammaEth+Y_F1*gammaB-Y_F1*gammaEth-gammaG
     +lu)/(Y_AE*alpha+2*Y_AE-1)/delt/2
         pd(1,5) = -1/Y_AE*miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))
         pd(4,3) = -(Y_AE*gammaEth+Y_F1*gammaB-Y_F1*gammaEth-gammaGlu)/g
     +ammaEth/(Y_AE*alpha+2*Y_AE-1)/delt/2
         pd(5,6) = 0
         pd(6,1) = -Y_F1*gammaEth*(Y_F2*gammaB-gammaGlu)/(Y_F1*gammaB*ga
     +mmaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaG
     +lu+gammaEth*gammaLac-gammaGlu*gammaLac)/delt
         pd(5,3) = Y_AE/(Y_AE*alpha+2*Y_AE-1)/delt/2
         pd(3,6) = -2*(Y_AE*alpha+2*Y_AE-1)/Y_F1*miuMax_F1*sp(1,j)/(k_sF
     +1+sp(1,j))
         pd(1,2) = 0
         pd(4,2) = 0
         pd(2,4) = 0
         pd(6,5) = 0
         pd(7,7) = -b_F2-1/delt
         pd(5,4) = 0
         pd(6,2) = 0
         pd(5,1) = 0
         pd(2,7) = 0
         pd(3,2) = 0
         pd(1,3) = 1/(Y_AE*alpha+2*Y_AE-1)/delt/2
         pd(3,7) = 0
         pd(6,4) = Y_F1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F
     +2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gam
     +maLac)*gammaEth*gammaLac/delt
         pd(2,2) = -1/Y_AE*miuMax_AE/(k_sAE+sp(2,j))*sp(5,j)+1/Y_AE*miuM
     +ax_AE*sp(2,j)/(k_sAE+sp(2,j))**2*sp(5,j)-1/delt
         pd(2,1) = 0
         pd(4,6) = 0
         pd(6,3) = Y_F1*(Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaLac)*ga
     +mmaEth/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*ga
     +mmaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)/(Y_A
     +E*alpha+2*Y_AE-1)/delt/2
         pd(7,6) = 0
         pd(1,6) = 0
         pd(5,2) = 0
         pd(5,7) = 0
         pd(2,3) = 0
      end
