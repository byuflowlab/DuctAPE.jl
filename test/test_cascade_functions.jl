@testset "Stagger+Inflow Cascade Types and Functions" begin
    @testset "Parse and Write File" begin
        info, stagger, Re, Mach, solidity, inflow, cl, cdrag = dt.parsecascadefile(
            "data/cas1111.dat", true
        )

        @test info == "Test Cascade File (made up data just for unit testing)"
        @test stagger == 0.0
        @test Re == 50000
        @test Mach == 0.0
        @test solidity == 0.25
        @test isapprox(inflow, [-pi; 0.0; pi])
        @test isapprox(cl, [0.0; 0.5; 0.0])
        @test isapprox(cdrag, [0.0; 0.05; 0.0])

        dt.writecascadefile(
            "data/cascadewritecheck.dat",
            info,
            stagger,
            inflow,
            Re,
            Mach,
            solidity,
            cl,
            cdrag,
            true,
        )

        info, stagger, Re, Mach, solidity, inflow, cl, cdrag = dt.parsecascadefile(
            "data/cascadewritecheck.dat", true
        )

        @test info == "Test Cascade File (made up data just for unit testing)"
        @test stagger == 0.0
        @test Re == 50000
        @test Mach == 0.0
        @test solidity == 0.25
        @test isapprox(inflow, [-pi; 0.0; pi])
        @test isapprox(cl, [0.0; 0.5; 0.0])
        @test isapprox(cdrag, [0.0; 0.05; 0.0])
    end

    @testset "Define Cascade" begin
        stagger = [0.0; 0.5; 1.0]
        inflow = [-pi; 0.0; pi]
        cl = repeat([0.0; 0.5; 0.0]; inner=(1, 3))
        cdrag = repeat([0.0; 0.05; 0.0]; inner=(1, 3))
        Re = 50000.0
        Mach = 0.0
        solidity = 0.25
        info = "Test Cascade File (made up data just for unit testing)"
        cas = dt.StaggerInflowCAS(
            "data/" .* ["cas1111.dat"; "cas2111.dat"; "cas3111.dat"]
        )
        st2 = dt.StaggerInflowCAS(stagger, inflow, cl, cdrag, info, Re, Mach, solidity)
        st3 = dt.StaggerInflowCAS(stagger, inflow, cl, cdrag)

        @test cas.stagger == stagger
        @test cas.inflow == inflow
        @test cas.cl == cl
        @test cas.cd == cdrag
        @test cas.info == info
        @test cas.Re == Re
        @test cas.Mach == Mach
        @test cas.solidity == solidity
        @test cas.stagger == st2.stagger
        @test cas.inflow == st2.inflow
        @test cas.cl == st2.cl
        @test cas.cd == st2.cd
        @test cas.info == st2.info
        @test cas.Re == st2.Re
        @test cas.Mach == st2.Mach
        @test cas.solidity == st2.solidity
        @test cas.stagger == st3.stagger
        @test cas.inflow == st3.inflow
        @test cas.cl == st3.cl
        @test cas.cd == st3.cd
        @test st3.info == "DuctTAPE written cascade"
        @test st3.Re == 0.0
        @test st3.Mach == 0.0
        @test st3.solidity == 0.0
    end

    @testset "Evaluate Cascade" begin
        cas = dt.StaggerInflowCAS(
            "data/" .* ["cas1111.dat"; "cas2111.dat"; "cas3111.dat"]
        )
        cl, cdrag = dt.caseval(cas, 0.5, 0.0;)

        @test cl == 0.5
        @test cdrag == 0.05

        cl, cdrag = dt.caseval(cas, 0.25, pi / 2;)

        @test isapprox(cl, 0.375)
        @test isapprox(cdrag, 0.0375)
    end

    @testset "Write Cascasde from Struct" begin
        cas = dt.StaggerInflowCAS(
            "data/" .* ["cas1111.dat"; "cas2111.dat"; "cas3111.dat"]
        )

        dt.writecascadefile(
            "data/" .* ["test1111.dat"; "test2111.dat"; "test3111.dat"],
            cas;
            radians=true,
        )

        cas2 = dt.StaggerInflowCAS(
            "data/" .* ["test1111.dat"; "test2111.dat"; "test3111.dat"]
        )

        @test cas.stagger == cas2.stagger
        @test cas.inflow == cas2.inflow
        @test cas.cl == cas2.cl
        @test cas.cd == cas2.cd
        @test cas.info == cas2.info
        @test cas.Re == cas2.Re
        @test cas.Mach == cas2.Mach
        @test cas.solidity == cas2.solidity
    end
end

@testset "Stagger+Inflow+Re+Mach Cascade Types and Functions" begin
    @testset "Define Cascade" begin
        stagger = [0.0; 0.5; 1.0]
        inflow = [-pi; 0.0; pi]
        cl = repeat([0.0; 0.5; 0.0]; inner=(1, 3, 3, 3))
        cdrag = repeat([0.0; 0.05; 0.0]; inner=(1, 3, 3, 3))
        Re = [50000.0; 100000; 150000]
        Mach = [0.0; 0.1; 0.2]
        solidity = 0.25
        filenames = Array{String}(undef, 3, 3, 3)
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    filenames[i, j, k] = "data/cas$i$j$(k)1.dat"
                end
            end
        end

        info = "Test Cascade File (made up data just for unit testing)"
        cas = dt.StaggerInflowReMachCAS(filenames)
        st2 = dt.StaggerInflowReMachCAS(
            stagger, inflow, Re, Mach, cl, cdrag, info, solidity
        )
        st3 = dt.StaggerInflowReMachCAS(stagger, inflow, Re, Mach, cl, cdrag)

        @test cas.stagger == stagger
        @test cas.inflow == inflow
        @test cas.cl == cl
        @test cas.cd == cdrag
        @test cas.info == info
        @test cas.Re == Re
        @test cas.Mach == Mach
        @test cas.solidity == solidity
        @test cas.stagger == st2.stagger
        @test cas.inflow == st2.inflow
        @test cas.cl == st2.cl
        @test cas.cd == st2.cd
        @test cas.info == st2.info
        @test cas.Re == st2.Re
        @test cas.Mach == st2.Mach
        @test cas.solidity == st2.solidity
        @test cas.stagger == st3.stagger
        @test cas.inflow == st3.inflow
        @test cas.cl == st3.cl
        @test cas.cd == st3.cd
        @test st3.info == "DuctTAPE written cascade"
        @test st3.Re == Re
        @test st3.Mach == Mach
        @test st3.solidity == 0.0
    end

    @testset "Evaluate Cascade" begin
        filenames = Array{String}(undef, 3, 3, 3, 3)
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    for ell in 1:3
                        filenames[i, j, k, ell] = "data/cas$i$j$k$ell.dat"
                    end
                end
            end
        end

        cas = dt.StaggerInflowReMachSolidityCAS(filenames)

        cl, cdrag = dt.caseval(cas, 0.5, pi, 100000, 0.1, 0.5)

        @test cl == 0.5
        @test cdrag == 0.05

        cl, cdrag = dt.caseval(cas, 0.25, pi / 2, 50000, 0.2, 0.75)

        @test isapprox(cl, 0.375)
        @test isapprox(cdrag, 0.0375)
    end

    @testset "Write Cascasde from Struct" begin end
end

@testset "Stagger+Inflow+Re+Mach+Solidity Cascade Types and Functions" begin
    @testset "Define Cascade" begin
        stagger = [0.0; 0.5; 1.0]
        inflow = [-pi; 0.0; pi]
        cl = repeat([0.0; 0.5; 0.0]; inner=(1, 3, 3, 3, 3))
        cdrag = repeat([0.0; 0.05; 0.0]; inner=(1, 3, 3, 3, 3))
        Re = [50000.0; 100000; 150000]
        Mach = [0.0; 0.1; 0.2]
        solidity = [0.25; 0.5; 0.75]
        info = "Test Cascade File (made up data just for unit testing)"

        # direct definitions
        st2 = dt.StaggerInflowReMachSolidityCAS(
            stagger, inflow, Re, Mach, solidity, cl, cdrag, info
        )
        st3 = dt.StaggerInflowReMachSolidityCAS(
            stagger, inflow, Re, Mach, solidity, cl, cdrag
        )

        # defined from input files.
        filenames = Array{String}(undef, 3, 3, 3, 3)
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    for ell in 1:3
                        filenames[i, j, k, ell] = "data/cas$i$j$k$ell.dat"
                    end
                end
            end
        end

        cas = dt.StaggerInflowReMachSolidityCAS(filenames)

        @test cas.stagger == stagger
        @test cas.inflow == inflow
        @test cas.cl == cl
        @test cas.cd == cdrag
        @test cas.info == info
        @test cas.Re == Re
        @test cas.Mach == Mach
        @test cas.solidity == solidity
        @test cas.stagger == st2.stagger
        @test cas.inflow == st2.inflow
        @test cas.cl == st2.cl
        @test cas.cd == st2.cd
        @test cas.info == st2.info
        @test cas.Re == st2.Re
        @test cas.Mach == st2.Mach
        @test cas.solidity == st2.solidity
        @test cas.stagger == st3.stagger
        @test cas.inflow == st3.inflow
        @test cas.cl == st3.cl
        @test cas.cd == st3.cd
        @test st3.info == "DuctTAPE written cascade"
        @test st3.Re == Re
        @test st3.Mach == Mach
        @test st3.solidity == solidity
    end

    @testset "Evaluate Cascade" begin
        filenames = Array{String}(undef, 3, 3, 3, 3)
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    for ell in 1:3
                        filenames[i, j, k, ell] = "data/cas$i$j$k$ell.dat"
                    end
                end
            end
        end

        cas = dt.StaggerInflowReMachSolidityCAS(filenames)

        cl, cdrag = dt.caseval(cas, 0.5, pi, 100000, 0.1, 0.5)

        @test cl == 0.5
        @test cdrag == 0.05

        cl, cdrag = dt.caseval(cas, 0.25, pi / 2, 50000, 0.2, 0.75)

        @test isapprox(cl, 0.375)
        @test isapprox(cdrag, 0.0375)
    end

    @testset "Write Cascasde from Struct" begin end
end
