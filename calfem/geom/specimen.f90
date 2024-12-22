SUBROUTINE SPECIMEN(Element,Nodes,ndoff,nelm,bc,nbc,loaddof,nldof,loaddof2,nldof2)

INTEGER                                      :: matOpen, matGetDir, matGetNextVariable,mp, pa
INTEGER                                      :: matGetNextVariableInfo, matGetVariable,mxCreateDoubleMatrix
INTEGER                                      :: i,j,ndoff,nelm,nelm_x,nelm_y,GT,n,elgrp,nldof2
INTEGER,             DIMENSION(nelm,9)       :: Element
DOUBLE PRECISION,    DIMENSION(nelm,9)       :: El_in
DOUBLE PRECISION,    DIMENSION(ndoff,2)      :: Nodes
DOUBLE PRECISION,    DIMENSION(nbc,2)        :: bc
INTEGER,             DIMENSION(nldof)        :: loaddof
DOUBLE PRECISION,    DIMENSION(nldof)        :: loaddof_in
INTEGER,             DIMENSION(nldof2)       :: loaddof2
DOUBLE PRECISION,    DIMENSION(nldof2)       :: loaddof2_in
INTEGER                                      :: mxGetM, mxGetN, matClose,ndir,  stat
!--------------------------------------------------------------
write(*,*)nbc


      mp = matOpen('geom.mat', 'r')

      pa=matGetVariable(mp, 'Node')      
      CALL mxCopyPtrToReal8(mxGetPr(pa), NODES, NDOFf*2)
      call mxDestroyArray(pa)
      
     
      
      pa=matGetVariable(mp, 'loaddof')
      CALL mxCopyPtrToReal8(mxGetPr(pa), loaddof_in, nldof)
      call mxDestroyArray(pa)
      
      pa=matGetVariable(mp, 'loaddof2')
      CALL mxCopyPtrToReal8(mxGetPr(pa), loaddof2_in, nldof2)
      call mxDestroyArray(pa)

      pa=matGetVariable(mp, 'element')
      CALL mxCopyPtrToReal8(mxGetPr(pa), el_in, NELM*9)
      call mxDestroyArray(pa)

      pa=matGetVariable(mp, 'bc')
      CALL mxCopyPtrToReal8(mxGetPr(pa), BC, 2*NBC)
      CALL mxDestroyArray(pa)


      stat = matClose(mp)

      element   =INT(el_in)
      loaddof   =INT(loaddof_in)
      loaddof2  =INT(loaddof2_in)
RETURN
END
