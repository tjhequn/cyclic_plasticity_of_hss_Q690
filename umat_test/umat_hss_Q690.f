C-----------------------------------------------------------------------
C     This umat is used for considering the cyclic softening behaviour 
C     of construction steel!
C                 written by HE Qun
C                   2020-05-16
C-----------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,cmname,ndi,nshr,ntens,
     3 nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)
C
      include 'aba_param.inc'
C
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),
     4 dfgrd0(3,3),dfgrd1(3,3)

C ----------------------------------------------------------------------
      integer, parameter :: num_alpha = 3, num_Y = 3

C material parameters declaration
      real*8 :: emod0, enu, sigma_y, strain_sh, sigma_sat, m_n, m_l, 
     1          m_beta, m_phi, phi_sat, m_alpha_list(num_alpha),  
     2          omega_list(num_alpha), m_Y_list(num_Y), Q_Y_list(num_Y),
     3          E_sat, delta_E, xi_E

C solution dependent variables declaration
      real*8 :: np_n(6), alpha_n(6), eqplas_n, Y_iso_n, phi_n, psi_n, 
     1          alpha_list_n(6,num_alpha), Y_iso_list_n(num_Y),
     2          alpha_mon_n(num_alpha),  
     3          np_n1(6), alpha_n1(6),eqplas_n1, Y_iso_n1, phi_n1, 
     4          psi_n1, alpha_list_n1(6,num_alpha),  
     5          Y_iso_list_n1(num_Y),
     6          alpha_mon_n1(num_alpha), g_n, g_n1, dg_n1

C other variables declaration
      real*8 :: matP1(6,6), matP2(6,6), matp3(6,6), s_trial(6), dlambda, 
     1          ss_trial(6), len_ss, xi_bar_n1(6), len_xi, dxi_dlm,  
     2          ss_n1(6), R_mon_n1, dR_mon_dp, dphi_dlm, dY_iso_dlm,    
     3          dlm_dstrain(6), dn_dstrain(6,6), 
     4          sum_alpha_1, sum_alpha_2, m_alpha_1(6),m_alpha_2(6), 
     4          temp_psi, m_alpha_mat(num_alpha,num_alpha), 
     5          sum_alpha_v(num_alpha)

C other parameters declaration     
      character*80 cmname
      parameter (max_iter=100, tolerance=1.0d-5)

C-----------------------------------------------------------------------
C material constants and parameters definition
C props(1 ):  Emod0, elastic modulus
C props(2 ):  enu, possion ratio
C props(3 ):  sigma_y, yield stress
C props(4 ):  strain_sh, length of yield plateau
C props(5 ):  sigma_sat, nonlinear hardening
C props(6 ):  m_n, nonlinear hardening modulus
C props(7 ):  m_l, linear hardening modulus
C props(8 ):  m_phi
C props(9):  phi_sat
C props(10:9+num_alpha):  m_alpha_list
C props(10+num_alpha:9+2*num_alpha):  omega_list
C props(10+2*num_alpha:9+2*num_alpha+num_Y):  m_Y_list
C props(10+2*num_alpha+num_Y:9+2*num_alpha+2*num_Y):  Q_Y_list
C props(10+2*num_alpha+2*num_Y):  E_sat
C props(11+2*num_alpha+2*num_Y):  xi_E
C-----------------------------------------------------------------------
C material constants and parameters assignment
C----------------------------------------------------------------------- 
      emod0        = props(1)
      enu          = props(2)
      sigma_y      = props(3)
      strain_sh    = props(4)
      sigma_sat    = props(5)
      m_n          = props(6)
      m_l          = props(7)
      m_phi        = props(8)
      phi_sat      = props(9)
      m_alpha_list = props(10:9+num_alpha)
      omega_list   = props(10+num_alpha:9+2*num_alpha)
      m_Y_list     = props(10+2*num_alpha:9+2*num_alpha+num_Y)
      Q_Y_list     = props(10+2*num_alpha+num_Y:9+2*num_alpha+2*num_Y)
      E_sat        = props(10+2*num_alpha+2*num_Y)
      xi_E         = props(11+2*num_alpha+2*num_Y)

      delta_E = 1.d0 - E_sat/Emod0
      
C ----------------------------------------------------------------------
C calculate shaer modulus and lamada
C ----------------------------------------------------------------------        
      eg = emod0/(2.d0*(1.d0+enu))
      elam = enu*emod0/((1.d0+enu)*(1.d0-2.d0*enu))

C-----------------------------------------------------------------------
C initialization of essential vector and matrix, such as elstic matrix 
C ddsdde, auxiliary matrix matP1, matP2, matP3
C matP1 is the matrix for inner product of stress
C matP2 is the matrix for the computation of deviatoric stress tensor
C matp3 is the matrix for the computation of deviatoric strain tensor
C-----------------------------------------------------------------------
      do k1=1, 6
          do k2=1, 6
              ddsdde(k2,k1) = 0.d0
              matP1(k2,k1)  = 0.d0
              matP2(k2,k1)  = 0.d0
              matP3(k2,k1)  = 0.d0
          end do
      end do

      do k1=1, 3
          do k2=1, 3
              ddsdde(k2, k1) = elam
              matP2(k2, k1) = -1.d0/3.d0
              matP3(k2, k1) = -1.d0/3.d0
          end do
          ddsdde(k1, k1) = 2.d0*eg+elam
          matP1(k1, k1) = 1.d0
          matP2(k1, k1) = 2.d0/3.d0
          matP3(k1, k1) = 2.d0/3.d0
      end do

      do k1=4, 6
          ddsdde(k1, k1) = eg
          matP1(k1, k1) = 2.d0
          matP2(k1, k1) = 1.d0
          matP3(k1, k1) = 1.d0/2.d0
      end do

      do k1=1,num_alpha
        sum_alpha_v(k1) = 1.d0
      end do

C-----------------------------------------------------------------------
C initialization of auxiliary constants
C-----------------------------------------------------------------------
      root32 = sqrt(3.d0/2.d0)
      root23 = 1.d0/root32

C-----------------------------------------------------------------------
C state variables definition
C statev(1)~statev(6): plastic flow direciton--np
C statev(7)~statev(12): backstress--alpha
C statev(13): equivalent plastic strain--eplas
C statev(14): isotropic hardening for cyclic load--Y_iso
C statev(15): backstress softening--phi
C statev(16): load mode indicator--psi
C statev(17)~statev(16+6*num_alpha): nonlinear backstress--alpha_list
C statev(17+6*num_alpha)~statev(16+6*num_alpha+num_Y): isotropic hardeni
C                                                      ng--Y_iso_list
C statev(17+6*num_alpha+num_Y)~statev(16+7*num_alpha+num_Y): alpha_n for 
C                                                    monotonic load mode
C statev(17+7*num_alpha+num_Y): E_mod_n1
C statev(18+7*num_alpha+num_Y): g_n1
C statev(19+7*num_alpha+num_Y): numer of iterations

C read the results of state variables from last step
C-----------------------------------------------------------------------
      np_n         = statev(1:6)
      alpha_n      = statev(7:12)
      eqplas_n     = statev(13)
      Y_iso_n      = statev(14)
      phi_n        = statev(15)
      psi_n        = statev(16)
      alpha_list_n = reshape(statev(17:16+6*num_alpha),[6,num_alpha])
      Y_iso_list_n = statev(17+6*num_alpha:16+6*num_alpha+num_Y)
      alpha_mon_n  = statev(17+6*num_alpha+num_Y:16+7*num_alpha+num_Y)
      g_n          = statev(18+7*num_alpha+num_Y)+1.d0

C-----------------------------------------------------------------------
C calculate the trial state
C-----------------------------------------------------------------------
C calculate the trial stress and its deviatoric part

      s_trial = stress + g_n*matmul(ddsdde,dstran)
      ss_trial = matmul(matP2,s_trial)

C calculate yield
      xi_bar_n1 = ss_trial - alpha_n
      len_xi = sqrt(dot_product(xi_bar_n1, matmul(matP1,xi_bar_n1)))
      f_n1_trial = root32*len_xi - sigma_y - Y_iso_n
      
C check the trial state
      if (f_n1_trial.le.tolerance) then

C update all the state variables
        stress = s_trial

        statev(1:6) = np_n
        statev(7:12) = alpha_n
        statev(13) = eqplas_n
        statev(14) = Y_iso_n
        statev(15) = phi_n
        statev(16) = psi_n
        statev(17:16+6*num_alpha)=reshape(alpha_list_n,[6*num_alpha])
        statev(17+6*num_alpha:16+6*num_alpha+num_Y) = Y_iso_list_n
        statev(17+6*num_alpha+num_Y:16+7*num_alpha+num_Y) = alpha_mon_n
        statev(17+7*num_alpha+num_Y) = g_n*emod0
        statev(18+7*num_alpha+num_Y) = g_n-1.d0
        statev(19+7*num_alpha+num_Y) = 0
        
        ddsdde = g_n*ddsdde
        

      else
C the trial state is not the real state, perform plastic correction
C initialize essential variables
        num_iter = 0
        dlambda = 0.0d0
        f_n1 = f_n1_trial

C-----------------------------------------------------------------------
C plastic correction
C-----------------------------------------------------------------------

        do while (abs(f_n1).gt.tolerance)
          num_iter = num_iter + 1

          if (num_iter.gt.max_iter) then

            pnewdt=0.25

          else

C--------------calculation of g_n1 and dg_n1 and dphi_dlm---------------
            eqplas_n1 = eqplas_n + dlambda
            g_n1 = 1.d0 - delta_E*(1.d0 - exp(-xi_E*eqplas_n1))
            dg_n1 = xi_E*(1.d0 - delta_E - g_n1)

C--------------calculation of g_n1 and dg_n1 and dphi_dlm---------------

C------------------calculation of phi_n1 and dphi_dlm-------------------
            m_alpha_1 =matmul(alpha_list_n,1.d0/(1.d0+
     1                      m_alpha_list*dlambda))
    
            Y_iso_list_n1 = (Y_iso_list_n+m_Y_list*Q_Y_list*dlambda)/
     1                      (1.d0 + m_Y_list*dlambda)
            Y_iso_n1 = sum(Y_iso_list_n1)
            dY_iso_dlm = sum(m_Y_list*(Q_Y_list-Y_iso_n)*dlambda/
     1                     (1.d0+m_Y_list*dlambda)**2)
    
            xi_bar_n1 = g_n1*ss_trial/g_n-beta_n-m_alpha_1
    
            len_xi=sqrt(dot_product(xi_bar_n1,matmul(matP1,xi_bar_n1)))
    
            np_n1 = xi_bar_n1/len_xi
    
            temp_psi = dot_product(np_n,matmul(matP1,np_n1))
    
            if ((psi_n.eq.0.d0).and.(temp_psi.lt.0.d0)) then
              psi_n1 = 1.d0
            else
              psi_n1 = psi_n
            end if
    
            if (psi_n1.eq.0.d0) then
    
              call monoHardening(sigma_y,strain_sh,sigma_sat,m_n,m_l,
     1                           eqplas_n1,R_mon_n1,dR_mon_dp)          
    
              temp_alpha_mono_n = sum(alpha_mon_n/(1.d0+
     1                                m_alpha_list*dlambda))
              dtemp_alpha_mono_n = sum(m_alpha_list*alpha_mon_n/
     1                                (1.d0+m_alpha_list*dlambda)**2)
    
              temp_F = sum(m_alpha_list*omega_list*dlambda/
     1                    (1.d0+m_alpha_list*dlambda))
              temp_G = R_mon_n1-temp_alpha_mono_n-sigma_y-Y_iso_n1
    
              dtemp_F = sum(m_alpha_list*omega_list/
     1                     (1.d0+m_alpha_list*dlambda)**2)
              dtemp_G = dR_mon_dp+dtemp_alpha_mono_n-dY_iso_dlm
    
              if (dlambda.eq.0.d0) then
                  phi_n1 = phi_n
                  dphi_dlm = 0.d0
              else
                  phi_n1 = temp_G/temp_F
                  dphi_dlm = (dtemp_G*temp_F-dtemp_F*temp_G)/temp_F**2
              end if
    
            else
    
              phi_n1 =(phi_n+m_phi*phi_sat*dlambda)/(1.d0+m_phi*dlambda)
              dphi_dlm = m_phi*(phi_sat-phi_n)*dlambda/(1.d0+
     1                  m_phi*dlambda)**2
            end if

C------------------calculation of phi_n1 and dphi_dlm-------------------
     
            sum_alpha_1 = sum(m_alpha_list*omega_list/(1.d0+
     1                        m_alpha_list*dlambda))
  
            sum_alpha_2 = sum(m_alpha_list*omega_list/(1.d0+
     1                        m_alpha_list*dlambda)**2)
  
            f_n1 = root32*len_xi-3.d0*g_n1*eg*dlambda-
     1             sum_alpha_1*phi_n1*dlambda-sigma_y-Y_iso_n1
  
            m_alpha_2 =matmul(alpha_list_n,m_alpha_list/(1.d0+
     1                        m_alpha_list*dlambda)**2)
  
            dalpha_dlm = sum_alpha_2*phi_n1+sum_alpha_1*dlambda*dphi_dlm
  
            dxi_dlm = dot_product(np_n1,matmul(matP1, 
     1                            dg_n1*ss_trial/g_n+m_alpha_2))
  
            df_n1 = root32*dxi_dlm-3.d0*g_n1*eg-3.d0*eg*dg_n1*dlambda-
     1              dalpha_dlm-dY_iso_dlm
  
            if (abs(f_n1).gt.tolerance) then
C update dlambda
              dlambda = dlambda - f_n1/df_n1

            else
C update related variables
              call diag(m_alpha_mat,1.d0/(1.d0+m_alpha_list*dlambda),
     1                  num_alpha)
              alpha_list_n1 = matmul(alpha_list_n+root23*phi_n1*dlambda*
     1                        matmul(reshape(np_n1,[6,1]),
     2                        reshape(omega_list*m_alpha_list,
     3                        [1,num_alpha])),m_alpha_mat)
                
              alpha_n1 = matmul(alpha_list_n1,sum_alpha_v)
              alpha_n1 = beta_n1+alpha_n1
  
              alpha_mon_n1 =(alpha_mon_n+m_alpha_list*omega_list*phi_n1*
     1                        dlambda)/(1.d0+m_alpha_list*dlambda)
  
              stress =g_n1*s_trial/g_n-sqrt(6.0d0)*g_n1*eg*dlambda*np_n1
  
              statev(1:6) = np_n1
              statev(7:12) = alpha_n1
              statev(13) = eqplas_n + dlambda
              statev(14) = Y_iso_n1
              statev(15) = phi_n1
              statev(16) = psi_n1
              statev(17:16+6*num_alpha)=reshape(alpha_list_n1,
     1                                             [6*num_alpha])
              statev(17+6*num_alpha:16+6*num_alpha+num_Y)=Y_iso_list_n1
              statev(17+6*num_alpha+num_Y:16+7*num_alpha+num_Y) = 
     1                                                     alpha_mon_n1
              statev(17+7*num_alpha+num_Y) = g_n1*emod0
              statev(18+7*num_alpha+num_Y) = g_n1-1.d0
              statev(19+7*num_alpha+num_Y) = num_iter


C compute tangential elastic modulus
              dlm_dstrain = -sqrt(6.0d0)*g_n1*eg*np_n1/df_n1
              dn_dstrain=2.d0*eg*g_n1*(matP3-matmul(
     1                   reshape(np_n1,[6,1]),reshape(np_n1,[1,6])))+
     2            matmul(reshape(dg_n1*ss_trial/g_n+m_alpha_2-
     3            dxi_dlm*np_n1,[6,1]),reshape(dlm_dstrain,[1,6]))
  
              dn_dstrain = dn_dstrain/len_xi
  
              ddsdde = ddsdde-sqrt(6.0d0)*eg*(dlambda*dn_dstrain+matmul(
     1                 reshape(np_n1,[6,1]),reshape(dlm_dstrain,[1,6])))
              ddsdde = g_n1*ddsdde + matmul(reshape(stress/g_n1,[6,1]), 
     2                             reshape(dg_n1*dlm_dstrain,[1,6]))

            end if

          end if 

        end do

      end if      
        
      return

      end
C------------------------------end of umat------------------------------


C-----------------------------------------------------------------------
C This subroutine is used to generate diagonal matrix
C-----------------------------------------------------------------------
      subroutine diag(A,b,N)

        real*8 :: A(N,N), b(N)

        do i=1,N
          do j=1,N
            A(i,j) = 0
          end do
          A(i,i) = b(i)
        end do

      return
      end
C------------------------------end of diag------------------------------


C-----------------------------------------------------------------------
C This subroutine is used for monotonic monoHardening
C-----------------------------------------------------------------------
      Subroutine monoHardening(sigma_y,strain_sh,sigma_sat,m_n,m_l,
     *                         eqplas,R_mon,dR_mon)

      real*8 :: sigma_y,strain_sh,sigma_sat,m_n,m_l,eqplas,R_mon,dR_mon

      if (eqplas.le.strain_sh) then
        R_mon = sigma_y
        dR_mon = 0.0d0
      else
        temp_peeq = eqplas-strain_sh
        R_mon=sigma_y+sigma_sat*(1.0d0-exp(-m_n*temp_peeq))+
     1                                                  m_l*temp_peeq
        dR_mon=sigma_sat*m_n*exp(-m_n*temp_peeq)+m_l
      end if
    
      return
      end
C--------------------------end of monoHardening-------------------------