#=
Utility functions for post-processing
=#

"""
    write_data(outs, filename; checkoutfileexists=false)

Writes NamedTuples, specifically for writing out the output of the `post_procces()` function.

# Arguments:
- `outs::NamedTuple` : Named tuple to write to file.
- `filename::String` : file name (including full desired path and file type) for file to write

# Keyword Arguments:
- `tuple_name::String` : desired variable name of written NamedTuple
- `checkoutfileexists::Bool=false` : boolean for whether to check if the outfile already exists and whether or not to overwrite it.

"""
function write_data(
    outs, filename="dfdc_outs.jl"; tuple_name="ductape_data", checkoutfileexists=false
)

    # check if file is present already
    if isfile(filename) && checkoutfileexists
        #ask if you want to overwrite it
        println(filename, " already exists. Overwrite?")
        println("n : no, exit")
        println("y or return : yes, continue")
        println("a : append the file")
        yn = readline(stdin)

        if any(contains.(yn, ["N", "n"]))
            #if no, exit
            return nothing
        elseif any(contains.(yn, ["A", "a"]))
            # if append, do so
            f = open(filename, "a")
            write(f, "$(tuple_name) = (;\n")
        else
            # otherwise overwrite
            f = open(filename, "w")
            write(f, "$(tuple_name) = (;\n")
        end
    else
        # overwrite if not checking
        f = open(filename, "w")
        write(f, "$(tuple_name) = (;\n")
    end

    # Write data to file
    # loop through keys and values
    for (k, p) in zip(keys(outs), outs)
        if typeof(p) <: String
            #need to put quotes on strings
            write(f, "\t$k = \"$(p)\",\n")
        else
            write(f, "\t$k = $p,\n")
        end
    end

    # close up
    write(f, ")")
    close(f)

    return nothing
end
