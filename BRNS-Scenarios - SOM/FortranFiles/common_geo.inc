        implicit real*8(a-h,o-z) 
        parameter (nsolid=0,ndiss=7,ncomp=7,nreac=7,nx=101)
        common/kinetics/miuMax_AE,k_sAE,miuMax_F1,k_sF1,miuMax_F2,k_sF2,
     +Y_AE,Y_F1,Y_F2,gammaGlu,gammaB,gammaEth,gammaLac,alpha,b_AE,b_F1,b
     +_F2
        real*8 miuMax_AE,k_sAE,miuMax_F1,k_sF1,miuMax_F2,k_sF2,Y_AE,Y_F1
     +,Y_F2,gammaGlu,gammaB,gammaEth,gammaLac,alpha,b_AE,b_F1,b_F2
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,endt
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,end
     +t
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic
        integer iq,iw,iDb,ipor,igrid,iarea,ic
