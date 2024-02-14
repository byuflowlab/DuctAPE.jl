      program test
        use integrate ! use integrate module containing LAMP and friends

C ----- Initialize variables and outputs
        REAL::X1,R1,X2,R2 
        REAL::XS,RS,XF,RF
        REAL::UG1,VG1,UG2,VG2
        REAL::US1,VS1,US2,VS2

C ------------------ C
C ----- TEST 2 ----- C
C ------------------ C
C ----- Assign Panel and Field Points
        X1 = 1.0
        R1 = 1.0
        X2 = 2.0
        R2 = 2.0
        XF = 2.0
        RF = 1.0
        XS = 1.5
        RS = 1.5

C ----- Open output Files
        open(unit=99,file="nominal_velocities1.jl")
        write(99,*) '# Nominal Velocities'
        write(99,*) 'p1 = [',X1,';',R1,']'
        write(99,*) 'p2 = [',X2,';',R2,']'
        write(99,*) 'pf = [',XF,';',RF,']'

        open(unit=88,file="self_velocities1.jl")
        write(88,*) '# Self-induced Velocities'
        write(88,*) 'p1 = [',X1,';',R1,']'
        write(88,*) 'p2 = [',X2,';',R2,']'
        write(88,*) 'ps = [',XS,';',RS,']'

C ----- Call Lamp and write outputs
        call LAMP(X1,R1,X2,R2,XF,RF,
     &            UG1,VG1,UG2,VG2,
     &            US1,VS1,US2,VS2)

        write(99,*) 'Vgammai = [',UG1,';',VG1,']'
        write(99,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(99,*) 'Vsigmai = [',US1,';',VS1,']'
        write(99,*) 'Vsigmaip1 = [',US2,';',VS2,']'
        
        open(unit=77,file="singular_components1.jl")
C ----- Call Lampc and write outputs
        call LAMPC(X1,R1,X2,R2,
     &             UG1,VG1,UG2,VG2,
     &             US1,VS1,US2,VS2)

        write(88,*) 'Vgammai = [',UG1,';',VG1,']'
        write(88,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(88,*) 'Vsigmai = [',US1,';',VS1,']'
        write(88,*) 'Vsigmaip1 = [',US2,';',VS2,']'
    
C ----- Close output files
        CLOSE(99)
        CLOSE(88)
        CLOSE(77)

C ------------------ C
C ----- TEST 2 ----- C
C ------------------ C
C ----- Assign Panel and Field Points
        X1 = 1.0
        R1 = 1.0
        X2 = 1.1
        R2 = 1.2
        XF = 3.0
        RF = 2.0
        XS = 1.05
        RS = 1.1

C ----- Open output Files
        open(unit=99,file="nominal_velocities2.jl")
        write(99,*) '# Nominal Velocities'
        write(99,*) 'p1 = [',X1,';',R1,']'
        write(99,*) 'p2 = [',X2,';',R2,']'
        write(99,*) 'pf = [',XF,';',RF,']'

        open(unit=88,file="self_velocities2.jl")
        write(88,*) '# Self-induced Velocities'
        write(88,*) 'p1 = [',X1,';',R1,']'
        write(88,*) 'p2 = [',X2,';',R2,']'
        write(88,*) 'ps = [',XS,';',RS,']'

C ----- Call Lamp and write outputs
        call LAMP(X1,R1,X2,R2,XF,RF,
     &            UG1,VG1,UG2,VG2,
     &            US1,VS1,US2,VS2)

        write(99,*) 'Vgammai = [',UG1,';',VG1,']'
        write(99,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(99,*) 'Vsigmai = [',US1,';',VS1,']'
        write(99,*) 'Vsigmaip1 = [',US2,';',VS2,']'
        
        open(unit=77,file="singular_components2.jl")
C ----- Call Lampc and write outputs
        call LAMPC(X1,R1,X2,R2,
     &             UG1,VG1,UG2,VG2,
     &             US1,VS1,US2,VS2)

        write(88,*) 'Vgammai = [',UG1,';',VG1,']'
        write(88,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(88,*) 'Vsigmai = [',US1,';',VS1,']'
        write(88,*) 'Vsigmaip1 = [',US2,';',VS2,']'
    
C ----- Close output files
        CLOSE(99)
        CLOSE(88)
        CLOSE(77)

C ------------------ C
C ----- TEST 3 ----- C
C ------------------ C
C ----- Assign Panel and Field Points
        X1 = 0.0
        R1 = 0.0
        X2 = 1.0
        R2 = 1.0
        XF = 3.0
        RF = 2.0
        XS = 0.5
        RS = 0.5

C ----- Open output Files
        open(unit=99,file="nominal_velocities3.jl")
        write(99,*) '# Nominal Velocities'
        write(99,*) 'p1 = [',X1,';',R1,']'
        write(99,*) 'p2 = [',X2,';',R2,']'
        write(99,*) 'pf = [',XF,';',RF,']'

        open(unit=88,file="self_velocities3.jl")
        write(88,*) '# Self-induced Velocities'
        write(88,*) 'p1 = [',X1,';',R1,']'
        write(88,*) 'p2 = [',X2,';',R2,']'
        write(88,*) 'ps = [',XS,';',RS,']'

C ----- Call Lamp and write outputs
        call LAMP(X1,R1,X2,R2,XF,RF,
     &            UG1,VG1,UG2,VG2,
     &            US1,VS1,US2,VS2)

        write(99,*) 'Vgammai = [',UG1,';',VG1,']'
        write(99,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(99,*) 'Vsigmai = [',US1,';',VS1,']'
        write(99,*) 'Vsigmaip1 = [',US2,';',VS2,']'
        
        open(unit=77,file="singular_components3.jl")
C ----- Call Lampc and write outputs
        call LAMPC(X1,R1,X2,R2,
     &             UG1,VG1,UG2,VG2,
     &             US1,VS1,US2,VS2)

        write(88,*) 'Vgammai = [',UG1,';',VG1,']'
        write(88,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(88,*) 'Vsigmai = [',US1,';',VS1,']'
        write(88,*) 'Vsigmaip1 = [',US2,';',VS2,']'
    
C ----- Close output files
        CLOSE(99)
        CLOSE(88)
        CLOSE(77)

C ------------------ C
C ----- TEST 4 ----- C
C ------------------ C
C ----- Assign Panel and Field Points
        X1 = 0.0
        R1 = 1.95943486e-17
        X2 = 0.000473018037
        R2 = 0.0122489417
        XS = 0.0002365090185
        RS = 0.0061244708500000095
        XF = 0.5
        RF = 0.5

C ----- Open output Files
        open(unit=99,file="nominal_velocities4.jl")
        write(99,*) '# Nominal Velocities'
        write(99,*) 'p1 = [',X1,';',R1,']'
        write(99,*) 'p2 = [',X2,';',R2,']'
        write(99,*) 'pf = [',XF,';',RF,']'

        open(unit=88,file="self_velocities4.jl")
        write(88,*) '# Self-induced Velocities'
        write(88,*) 'p1 = [',X1,';',R1,']'
        write(88,*) 'p2 = [',X2,';',R2,']'
        write(88,*) 'ps = [',XS,';',RS,']'

C ----- Call Lamp and write outputs
        call LAMP(X1,R1,X2,R2,XF,RF,
     &            UG1,VG1,UG2,VG2,
     &            US1,VS1,US2,VS2)

        write(99,*) 'Vgammai = [',UG1,';',VG1,']'
        write(99,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(99,*) 'Vsigmai = [',US1,';',VS1,']'
        write(99,*) 'Vsigmaip1 = [',US2,';',VS2,']'

        open(unit=77,file="singular_components4.jl")
C ----- Call Lampc and write outputs
        call LAMPC(X1,R1,X2,R2,
     &             UG1,VG1,UG2,VG2,
     &             US1,VS1,US2,VS2)

        write(88,*) 'Vgammai = [',UG1,';',VG1,']'
        write(88,*) 'Vgammaip1 = [',UG2,';',VG2,']'
        write(88,*) 'Vsigmai = [',US1,';',VS1,']'
        write(88,*) 'Vsigmaip1 = [',US2,';',VS2,']'

C ----- Close output files
        CLOSE(99)
        CLOSE(88)
        CLOSE(77)
       end program test