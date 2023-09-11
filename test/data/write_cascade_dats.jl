
stagger = [0.0; 0.5; 1.0]
Re = [50000; 100000; 150000]
Ma = [0.0; 0.1; 0.2]
solidity = [0.25; 0.5; 0.75]

for (i, stag) in enumerate(stagger)
    for (j, r) in enumerate(Re)
        for (k, m) in enumerate(Ma)
            for (ell, sol) in enumerate(solidity)
                f = open("test/data/cas$i$j$k$ell.dat", "w")
                write(f, "Test Cascade File (made up data just for unit testing)\n")
                write(f, "$stag\n")
                write(f, "$r\n")
                write(f, "$m\n")
                write(f, "$sol\n")
                write(f, "-3.1415926535897931 0.0 0.0\n")
                write(f, "0.0 0.5 0.05\n")
                write(f, "3.1415926535897931 0.0 0.0")
                close(f)
            end
        end
    end
end
