# These functions are R implementations of Python implementations of c++
# code required to correctly transform Siemens MRS geometry parameters to NIfTI
# affine format. Credit to Will Clarke (https://github.com/wexeee/) for figuring
# this out and making the Python code available at:
# https://github.com/wexeee/spec2nii/blob/master/spec2nii/GSL/gslfunctions.py

# R implementation of:
# IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLAlmEqual.cpp
fGSLAlmEqual <- function(dArg1, dArg2) {
  dTmp <- dArg1 - dArg2
  bOut <- (dTmp >= -1.0e-6)  &  (dTmp <= 1.0e-6)
  return(bOut)
}

# R implementation of:
# IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLClassOri.cpp
#
# * This function determines whether a supplied normal vector
# *	              describes a sagittal, coronal or transverse slice.
# *		Result:
# *		  CASE = 0: Sagittal
# *		  CASE = 1: Coronal
# *		  CASE = 2: Transverse
# *     This was created with the following table:
# *			(Note: ~ means 'about equal')
# *
# *		  |sag|~|cor| and |sag|~|tra|  -->  tra
# *
# *		  |sag|~|cor| and |sag|<|tra|  -->  tra
# *		  |sag|~|cor| and |sag|>|tra|  -->  cor
# *
# *		  |sag|~|tra| and |sag|<|cor|  -->  cor
# *		  |sag|~|tra| and |sag|>|cor|  -->  tra
# *
# *		  |cor|~|tra| and |cor|<|sag|  -->  sag
# *		  |cor|~|tra| and |cor|>|sag|  -->  tra
# *
# *		  |sag|>|cor| and |sag|>|tra|  -->  sag
# *		  |sag|>|cor| and |sag|<|tra|  -->  tra
# *		  |sag|<|cor| and |cor|>|tra|  -->  cor
# *		  |sag|<|cor| and |cor|<|tra|  -->  tra
# *
# *		  |sag|>|tra| and |sag|<|cor|  -->  cor
# *		  |sag|>|tra| and |sag|>|cor|  -->  sag
# *		  |sag|<|tra| and |tra|<|cor|  -->  cor
# *		  |sag|<|tra| and |tra|>|cor|  -->  tra
# *
# *		  |cor|>|tra| and |cor|<|sag|  -->  sag
# *		  |cor|>|tra| and |cor|>|sag|  -->  cor
# *		  |cor|<|tra| and |tra|<|sag|  -->  sag
# *		  |cor|<|tra| and |tra|>|sag|  -->  tra
# In:
# double  dSagComp,	   /* IMP: Sagittal component of normal vector          */
# double  dCorComp,	   /* IMP: Coronal component of normal vector           */
# double  dTraComp,	   /* IMP: Transverse component of normal vector        */
# Out:
# long *  plCase       /* EXP: Case (0=Sagittal, 1=Coronal or 2=Transverse) */
fGSLClassOri <- function(dSagComp, dCorComp, dTraComp, DEBUG = FALSE) {
  if (DEBUG) {
    print(sprintf('Normal vector = {dSagComp : %10.7f} {dCorComp : %10.7f} {dTraComp : %10.7f}.',
                     dSagComp, dCorComp, dTraComp))
  }
  
  #Compute some temporary values
  dAbsSagComp     <- abs(dSagComp)
  dAbsCorComp     <- abs(dCorComp)
  dAbsTraComp     <- abs(dTraComp)
  
  bAlmEqualSagCor <- fGSLAlmEqual(dAbsSagComp, dAbsCorComp)
  bAlmEqualSagTra <- fGSLAlmEqual(dAbsSagComp, dAbsTraComp)
  bAlmEqualCorTra <- fGSLAlmEqual(dAbsCorComp, dAbsTraComp)
  
  # Check all values to determine the slice orientation (sag, cor, tra)
  if ((bAlmEqualSagCor              &  bAlmEqualSagTra)             |
      (bAlmEqualSagCor              &  (dAbsSagComp < dAbsTraComp)) |
      (bAlmEqualSagTra              &  (dAbsSagComp > dAbsCorComp)) |
      (bAlmEqualCorTra              &  (dAbsCorComp > dAbsSagComp)) |
      ((dAbsSagComp > dAbsCorComp)  &  (dAbsSagComp < dAbsTraComp)) |
      ((dAbsSagComp < dAbsCorComp)  &  (dAbsCorComp < dAbsTraComp)) |
      ((dAbsSagComp < dAbsTraComp)  &  (dAbsTraComp > dAbsCorComp)) |
      ((dAbsCorComp < dAbsTraComp)  &  (dAbsTraComp > dAbsSagComp))) {
    
    #Mainly transverse...
    if (DEBUG) print('Mainly transverse.')        
    plcase <- 2 #CASE = 2: Transverse
  } else if ((bAlmEqualSagCor             &  (dAbsSagComp > dAbsTraComp)) |
             (bAlmEqualSagTra             &  (dAbsSagComp < dAbsCorComp)) |
            ((dAbsSagComp < dAbsCorComp)  &  (dAbsCorComp > dAbsTraComp)) |
            ((dAbsSagComp > dAbsTraComp)  &  (dAbsSagComp < dAbsCorComp)) |
            ((dAbsSagComp < dAbsTraComp)  &  (dAbsTraComp < dAbsCorComp))) {
    
    # Mainly coronal...
    if (DEBUG) print('Mainly coronal.')        
    plcase <- 1 # CASE = 1: Coronal
  } else if ((bAlmEqualCorTra             &  (dAbsCorComp < dAbsSagComp)) |
            ((dAbsSagComp > dAbsCorComp)  &  (dAbsSagComp > dAbsTraComp)) |
            ((dAbsCorComp > dAbsTraComp)  & (dAbsCorComp < dAbsSagComp))  |
            ((dAbsCorComp < dAbsTraComp)  &  (dAbsTraComp < dAbsSagComp))) { 
    
    # Mainly sagittal...
    if (DEBUG) print('Mainly Sagittal.')
    plcase <- 0 # CASE = 0: Sagittal
  } else {
    # Invalid slice orientation...
    stop('Error: Invalid slice orientation')    
  }
  return(plcase)
}

# R implementation of: 
# IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLCalcPRS.cpp
# * Description : calculates the two vectors of phase encoding and and
# *               readout direction.                                 
# *               The calculation depends on the slice orientation (the slice
# *               normal vector) and on the angle of rotation around the s axis.
# *		            Every vector (Gp, Gr and Gs) has the three components sag, cor
# *               and tra, the description is patient oriented. All three
# *               vectors have a length of 1.0. The biggest component of Gs must
# *               have a positive sign.
# *
# *		Formulas for the rotation around the vektor s:
# *		    (a = cos (dPhi), b = sin (dPhi))
# *
# *		    new             old               rotation     base
# *		    vector  =   coordinate system   *  matrix    * vector
# *
# *		    (P_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (1)
# *		    (P_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
# *		    (P_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (0)
# *
# *		    (R_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
# *		    (R_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (1)
# *		    (R_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (0)
# *
# *		    (S_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
# *		    (S_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
# *		    (S_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (1)
# *
# *		    This multiplied:
# *
# *		    (P_sag)   (a * P_sag - b * R_sag)
# *		    (P_cor) = (a * P_cor - b * R_cor)
# *		    (P_tra)   (a * P_tra - b * R_tra)
# *
# *		    (R_sag)   (b * P_sag + a * R_sag)
# *		    (R_cor) = (b * P_cor + a * R_cor)
# *		    (R_tra)   (b * P_tra + a * R_tra)
# *
# *		    (S_sag)   (S_sag)
# *		    (S_cor) = (S_cor)	well!
# *		    (S_tra)   (P_tra)                
# %
# Input:
#   dGs: The GS vector (= slice normal vector)
#   dPhi: The rotational angle around Gs
# Output: 
#   dGp: The Gp vector
#   dGr: The Gr vector
# From libGSL.h
fGSLCalcPRS <- function(dGs, dPhi, DEBUG = FALSE) {
  # patient axes
  SAGITTAL   <- 0
  CORONAL    <- 1
  TRANSVERSE <- 2
  
  # Start of function
  lOrientation <- 0  # Orientation (SAGITTAL, CORONAL or TRANSVERSE)
  lOrientation <- lOrientation + fGSLClassOri(dGs[SAGITTAL   + 1],
                                              dGs[CORONAL    + 1],
                                              dGs[TRANSVERSE + 1], DEBUG)
  dGp <- rep(0, 3)
  
  if (lOrientation == TRANSVERSE) {
    dGp[1] <-  0.0
    dGp[2] <-  dGs[3] * sqrt(1. / (dGs[2] * dGs[2] + dGs[3] * dGs[3]))
    dGp[3] <- -dGs[2] * sqrt(1. / (dGs[2] * dGs[2] + dGs[3] * dGs[3]))
  } else if (lOrientation == CORONAL) {
    dGp[1] <-  dGs[2] * sqrt(1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
    dGp[2] <- -dGs[1] * sqrt(1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
    dGp[3] <-  0.0
  } else if (lOrientation == SAGITTAL) {
    dGp[1] <- -dGs[2] * sqrt(1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
    dGp[2] <-  dGs[1] * sqrt(1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
    dGp[3] <-  0.0
  } else {
    stop('Invalid slice orientation returned from fGSLClassOri')
  }
  
  # Calculate GR = GS x GP 
  dGr <- rep(0, 3)
  dGr[1] <- dGs[2] * dGp[3] - dGs[3] * dGp[2]
  dGr[2] <- dGs[3] * dGp[1] - dGs[1] * dGp[3]
  dGr[3] <- dGs[1] * dGp[2] - dGs[2] * dGp[1]
  
  if (DEBUG) {
    print('Before rotation around S:')
    print(sprintf('GP = %10.7f %10.7f %10.7f', dGp[1], dGp[2], dGp[3]))
    print(sprintf('GR = %10.7f %10.7f %10.7f', dGr[1], dGr[2], dGr[3]))
    print(sprintf('GS = %10.7f %10.7f %10.7f', dGs[1], dGs[2], dGs[3]))
  }
  
  # Rotation
  if (dPhi != 0.0) {
    #Rotate around the S axis                                             
    if (DEBUG) {
      tmp <- dPhi * 180.0 / pi
      print(sprintf('PHI = %10.7f', tmp))
    }
  }
  
  dGp[1] <- cos(dPhi) * dGp[1] - sin(dPhi) * dGr[1]
  dGp[2] <- cos(dPhi) * dGp[2] - sin(dPhi) * dGr[2]
  dGp[3] <- cos(dPhi) * dGp[3] - sin(dPhi) * dGr[3]
  
  # Calculate new GR = GS x GP                                           
  dGr[1] <- dGs[2] * dGp[3] - dGs[3] * dGp[2]
  dGr[2] <- dGs[3] * dGp[1] - dGs[1] * dGp[3]
  dGr[3] <- dGs[1] * dGp[2] - dGs[2] * dGp[1]
  
  if (DEBUG) {
    print('After the Rotation around S:')
    print(sprintf('GP = %10.7f %10.7f %10.7f', dGp[1], dGp[2], dGp[3]))
    print(sprintf('GR = %10.7f %10.7f %10.7f', dGr[1], dGr[2], dGr[3]))
    print(sprintf('GS = %10.7f %10.7f %10.7f', dGs[1], dGs[2], dGs[3]))
  }
  
  return(list(dGp = dGp, dGr = dGr))
}