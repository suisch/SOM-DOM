c      
c     SUBROUTINE residual
c      
      subroutine residual(funcs,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension funcs(ncomp)
        call switches(j)
          funcs(1) = -1/Y_AE*miuMax_AE*sp(1,j)/(k_sAE+sp(1,j))*sp(5,j)+1
     +/(Y_AE*alpha+2*Y_AE-1)*sp(3,j)/delt/2-1/(Y_AE*alpha+2*Y_AE-1)*spol
     +d(3,j)/delt/2
          funcs(2) = -1/Y_AE*miuMax_AE*sp(2,j)/(k_sAE+sp(2,j))*sp(5,j)-s
     +p(2,j)/delt+spold(2,j)/delt
          funcs(3) = -2*(Y_AE*alpha+2*Y_AE-1)/Y_F1*miuMax_F1*sp(1,j)/(k_
     +sF1+sp(1,j))*sp(6,j)-(2*gammaEth*(Y_AE*Y_F2*alpha*gammaB+2*Y_AE*Y_
     +F2*gammaB-Y_AE*alpha*gammaGlu-2*Y_AE*gammaGlu-Y_F2*gammaB+gammaGlu
     +)/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEt
     +h-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*sp(1,j)-(
     +Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaLac)*gammaEth/(Y_F1*gammaB
     +*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gam
     +maGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*sp(3,j)-2*(Y_AE*alpha+2
     +*Y_AE-1)/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*
     +gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*ga
     +mmaEth*gammaLac*sp(4,j))/delt+(2*gammaEth*(Y_AE*Y_F2*alpha*gammaB+
     +2*Y_AE*Y_F2*gammaB-Y_AE*alpha*gammaGlu-2*Y_AE*gammaGlu-Y_F2*gammaB
     ++gammaGlu)/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gamma
     +B*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*
     +spold(1,j)-(Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaLac)*gammaEth/
     +(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-
     +gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*spold(3,j)-
     +2*(Y_AE*alpha+2*Y_AE-1)/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaL
     +ac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaG
     +lu*gammaLac)*gammaEth*gammaLac*spold(4,j))/delt
          funcs(4) = -(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*
     +gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gamma
     +Lac)/gammaEth/Y_F2/gammaLac*miuMax_F2*sp(1,j)/(k_sF2+sp(1,j))*sp(7
     +,j)-((Y_F1*gammaB-Y_F1*gammaEth+gammaEth-gammaGlu)/gammaEth*sp(1,j
     +)+(Y_AE*gammaEth+Y_F1*gammaB-Y_F1*gammaEth-gammaGlu)/gammaEth/(Y_A
     +E*alpha+2*Y_AE-1)*sp(3,j)/2+sp(4,j))/delt+((Y_F1*gammaB-Y_F1*gamma
     +Eth+gammaEth-gammaGlu)/gammaEth*spold(1,j)+(Y_AE*gammaEth+Y_F1*gam
     +maB-Y_F1*gammaEth-gammaGlu)/gammaEth/(Y_AE*alpha+2*Y_AE-1)*spold(3
     +,j)/2+spold(4,j))/delt
          funcs(5) = -b_AE*sp(5,j)-(-Y_AE/(Y_AE*alpha+2*Y_AE-1)*sp(3,j)/
     +2+sp(5,j))/delt+(-Y_AE/(Y_AE*alpha+2*Y_AE-1)*spold(3,j)/2+spold(5,
     +j))/delt
          funcs(6) = -b_F1*sp(6,j)-(Y_F1*gammaEth*(Y_F2*gammaB-gammaGlu)
     +/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth
     +-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)*sp(1,j)-Y_
     +F1*(Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaLac)*gammaEth/(Y_F1*ga
     +mmaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth
     +*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac)/(Y_AE*alpha+2*Y_AE-
     +1)*sp(3,j)/2-Y_F1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F
     +2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gam
     +maLac)*gammaEth*gammaLac*sp(4,j)+sp(6,j))/delt+(Y_F1*gammaEth*(Y_F
     +2*gammaB-gammaGlu)/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_
     +F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*ga
     +mmaLac)*spold(1,j)-Y_F1*(Y_AE*gammaLac-Y_F2*gammaB+gammaGlu-gammaL
     +ac)*gammaEth/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*gam
     +maB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLac
     +)/(Y_AE*alpha+2*Y_AE-1)*spold(3,j)/2-Y_F1/(Y_F1*gammaB*gammaLac-Y_
     +F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaE
     +th*gammaLac-gammaGlu*gammaLac)*gammaEth*gammaLac*spold(4,j)+spold(
     +6,j))/delt
          funcs(7) = -b_F2*sp(7,j)-(1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEt
     +h*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLa
     +c-gammaGlu*gammaLac)*Y_F2*gammaLac*(Y_F1*gammaB-Y_F1*gammaEth+gamm
     +aEth-gammaGlu)*sp(1,j)+1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gamma
     +Lac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gamma
     +Glu*gammaLac)*Y_F2*gammaLac*(Y_AE*gammaEth+Y_F1*gammaB-Y_F1*gammaE
     +th-gammaGlu)/(Y_AE*alpha+2*Y_AE-1)*sp(3,j)/2+1/(Y_F1*gammaB*gammaL
     +ac-Y_F1*gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+g
     +ammaEth*gammaLac-gammaGlu*gammaLac)*gammaEth*Y_F2*gammaLac*sp(4,j)
     ++sp(7,j))/delt+(1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F
     +2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gam
     +maLac)*Y_F2*gammaLac*(Y_F1*gammaB-Y_F1*gammaEth+gammaEth-gammaGlu)
     +*spold(1,j)+1/(Y_F1*gammaB*gammaLac-Y_F1*gammaEth*gammaLac+Y_F2*ga
     +mmaB*gammaEth-gammaEth*gammaGlu+gammaEth*gammaLac-gammaGlu*gammaLa
     +c)*Y_F2*gammaLac*(Y_AE*gammaEth+Y_F1*gammaB-Y_F1*gammaEth-gammaGlu
     +)/(Y_AE*alpha+2*Y_AE-1)*spold(3,j)/2+1/(Y_F1*gammaB*gammaLac-Y_F1*
     +gammaEth*gammaLac+Y_F2*gammaB*gammaEth-gammaEth*gammaGlu+gammaEth*
     +gammaLac-gammaGlu*gammaLac)*gammaEth*Y_F2*gammaLac*spold(4,j)+spol
     +d(7,j))/delt
      end
