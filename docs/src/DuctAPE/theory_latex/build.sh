FNAME=theory

# OPTION `build.sh genall`: Generate all auxiliary files
if [ ${1:-not} = genall ]
then

    # Remove all auxiliary files
    sh clean.sh

    # # Remove biber cache because it gets corrupted now and then.
    # rm -rf $(biber --cache)

    # Pre-compilation
    pdflatex ${FNAME}.tex

    # Generate index files
    makeindex ${FNAME}.nlo -s nomencl.ist -o ${FNAME}.nls
    # makeindex main.idx -s main.mst

    # Generate .blg .bbl bibliography files
    biber ${FNAME}

fi

# Pre-compile again to add bibliography to TOC
pdflatex ${FNAME}.tex

# Final compilation
pdflatex ${FNAME}.tex

# Remove auxiliary files
sh clean.sh
