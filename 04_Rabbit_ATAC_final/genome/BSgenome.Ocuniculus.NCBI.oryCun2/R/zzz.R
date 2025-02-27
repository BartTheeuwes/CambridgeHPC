###
###

.pkgname <- "BSgenome.Ocuniculus.NCBI.oryCun2"

.seqnames <- NULL

.circ_seqs <- character(0)

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Oryctolagus cuniculus",
        common_name="Rabbit",
        genome="OryCun2.0",
        provider="NCBI",
        release_date="Oct. 2009",
        source_url="https://www.ncbi.nlm.nih.gov/assembly/GCF_000003625.3",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Ocuniculus"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

