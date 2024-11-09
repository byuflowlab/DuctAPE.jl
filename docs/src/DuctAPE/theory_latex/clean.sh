FNAME=*

# OPTION `clean.sh keepnombib`: Keep auxiliary nomenclature and bibliography files
if [ ${1:-all} = keepaux ]; then
    EXTS="synctex.gz out log brf fls fdb_latexmk run.xml dvi nlg bcf nlo ilg"

# OPTION `clean.sh keepbib`: Keep auxiliary bibliography files
elif [ ${1:-all} = keepbib ]; then
    EXTS="synctex.gz out log brf fls fdb_latexmk run.xml dvi nlg bcf nlo ilg nls aux toc lof lot"

# DEFAULT `clean.sh`: Delete all auxiliary files
else
    EXTS="synctex.gz out log brf fls fdb_latexmk run.xml dvi nlg bcf nlo ilg nls aux toc lof lot blg bbl tdo ind idx mst auxlock"
fi

# Remove files
for EXT in ${EXTS}
    do
        rm -f ${FNAME}.${EXT}
done
