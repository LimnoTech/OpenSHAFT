        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 14 11:19:47 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETLOOKUPVALUE__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE GETLOOKUPVALUE(LOOKUPIN,LOOKUPOUT,     &
     &TABLEID,XCOL,YCOL)
              REAL(KIND=8), INTENT(IN) :: LOOKUPIN
              REAL(KIND=8), INTENT(INOUT) :: LOOKUPOUT
              INTEGER(KIND=4), INTENT(IN) :: TABLEID
              INTEGER(KIND=4), INTENT(IN) :: XCOL
              INTEGER(KIND=4), INTENT(IN) :: YCOL
            END SUBROUTINE GETLOOKUPVALUE
          END INTERFACE 
        END MODULE GETLOOKUPVALUE__genmod
