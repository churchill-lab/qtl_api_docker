#
# ASSUMING that the Rdata file has already been loaded
#


# 
# Load the libraries will we need
#

library(data.table)
library(dplyr)
library(intermediate)
library(jsonlite)
library(qtl2convert)
library(qtl2geno)
library(qtl2plot)
library(qtl2scan)
library(RSQLite)
library(pryr)

#
# Define some global variables
#

db.file <- "/api/data/ccfoundersnps.sqlite"

# Calculate the probs and create a map
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")


# we no longer are using genoprobs so free up some memory
genoprobs_env <- where("genoprobs")
rm(genoprobs, envir=genoprobs_env)
gc()

#
# fix the value passed in via the web
#  
fix_boolean <- function(b) {
    if (is.logical(b)) {
        return (b)
    } else if (is.character(b)) {
        return (toupper(b) %in% c("T", "TRUE", "YES", "Y", "1"))
    } else if (is.numeric(b)) {
        return (b == 1)
    }
    
    return (FALSE)
}

#
# fix numeric values with a default if it cannot 
#
nvl_integer <- function(n, default) {
    result <- tryCatch(
        {
            n <- as.numeric(n)
            if ((n %% 1) == 0) {
                return (n)
            } else {
                return (default)
            }
        },
        error=function(cond) {
            #message(cond)
            return (default)
        },
        warning=function(cond) {
            #message(cond)
            return (default)
        },
        finally={
            # nothing
        }
    )
    return (result)
}
#
# The following is used for debugging to show where the requests are coming from
#

#* @filter logger
function(req){
    #print(paste0("HTTP_ACCEPT: ", req$HTTP_ACCEPT))
    #print(paste0("HTTP_ACCEPT_ENCODING: ", req$HTTP_ACCEPT_ENCODING))
    #print(paste0("HTTP_CONNECTION: ", req$HTTP_CONNECTION))
    #print(paste0("HTTP_HOST: ", req$HTTP_HOST))
    #print(paste0("HTTP_USER_AGENT: ", req$HTTP_USER_AGENT))
    #print(paste0("PATH_INFO: ", req$PATH_INFO))
    #print(paste0("QUERY_STRING: ", req$QUERY_STRING))
    #print(paste0("REMOTE_ADDR: ", req$REMOTE_ADDR))
    #print(paste0("REMOTE_PORT: ", req$REMOTE_PORT))
    #print(paste0("REQUEST_METHOD: ", req$REQUEST_METHOD))
    #print(paste0("SCRIPT_NAME: ", req$SCRIPT_NAME))
    #print(paste0("SERVER_NAME: ", req$SERVER_NAME))
    #print(paste0("SERVER_PORT: ", req$SERVER_PORT))
    
    print(paste0(date(), " - ",
                 req$REMOTE_ADDR, " - ",
                 req$REQUEST_METHOD, " ",
                 req$PATH_INFO, req$QUERY_STRING))
    forward()
}

#
# show the time it took for the request 
#
track_time <- function(req, elapsed) {
    print(sprintf("URL: [%s%s] TIME: [%.3f]", req$PATH_INFO, req$QUERY_STRING, elapsed))
}    

#
# Generate an error response 
#
set_error <- function(res, code=400, message="Error") {
    res$status <- code
    res$body <- toJSON(list(message=message), auto_unbox = TRUE)
    res
}

#' Get the system info
#'
#' @param req the request object
#' @param res the response object
#'
#' @return JSON of the system information
#'
#* @get /info
http_info <- function(req, res) {
    # start the clock
    ptm <- proc.time()

    to_return <- as.list(Sys.info())
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])

    return (to_return)
}

#' Get the options
#'
#' @param req the request object
#' @param res the response object
#'
#' @return JSON object of the options
#'
#' Example:
#'    {"data_levels":["mrna"],
#'     "covar_factors":[{"column_name":"Diet","display_name":"Diet"},
#'                      {"column_name":"Sex","display_name":"Sex"}
#'                     ]
#'    }
#'
#* @get /options
http_options <- function(req, res) {
    # start the clock
    ptm <- proc.time()

    data_levels <- ''
    if (exists('annot.mrna')) {
        data_levels <- c('mrna') 
    }
    if (exists('annot.protein')) {
        data_levels <- c(data_levels, 'protein') 
    }

    to_return <- list(data_levels=data_levels, covar_factors=covar_factors)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])

    return (to_return)
}

#' Perform the LOD scan
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl genen identifier
#' @param regress_local TRUE to regress on local genotype, FALSE to not
#' @param ncores number of cores to use (0=ALL)
#'
#' @return JSON data of LOD scan
#'
#' Each element in the array is:
#'
#' id: string (Ensembl ID)
#' chr: string
#' pos: integer
#' pheno1: numeric
#'
#' Example:
#'     {"data":[{"1_40055","1",40055,2.4787},
#'              {"1_87283","1",87283,2.4787},
#'              ...
#'              {"1_134512","1",134512,:2.4787}
#'             ]}
#'
#* @get /lodscan/mrna
#* @post /lodscan/mrna
http_lodscan_mrna <- function(req, res, id, regress_local=FALSE, ncores=0) {
    # start the clock
    ptm <- proc.time()
    
    # get the Ensembl gene index 
    idx <- which(annot.mrna$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }
    
    # make sure ncores is appropriate  
    num_cores = nvl_integer(ncores, 0)
    
    # set covars
    addcovar <- covar[,-1]
    
    # to regress local genotype, add neareast marker to covariates
    if (fix_boolean(regress_local)) {
        #addcovar <- cbind(addcovar, genoprobs[,-1,annot.mrna$nearest_snp[idx]])
        mkr_tmp = as.character(snps[annot.mrna$nearest_snp[idx],1])
        chr_tmp = as.character(snps[annot.mrna$nearest_snp[idx],2])
        addcovar <- cbind(addcovar, probs[,chr_tmp][[1]][,-1,mkr_tmp])
    }
 
    # perform the scan using QTL2
    temp <- (scan1(genoprobs=probs, 
                   kinship=Glist, 
                   pheno=expr.mrna[,idx], 
                   addcovar=addcovar, 
                   cores=num_cores, 
                   reml=TRUE))
    
    # construct a 2 dimensional array of data
    temp_data <- data.table(id=snps$marker, chr=snps$chr, pos=snps$pos, temp)
    # setting colnames to NULL removes the names in the JSON and return an array
    colnames(temp_data) <- NULL
    to_return <- list(data=temp_data)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])

    return (to_return)
}

#' Perform the LOD scan
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl protein identifier
#' @param regress_local TRUE to regress on local genotype, FALSE to not
#' @param ncores number of cores to use (0=ALL)
#'
#' @return JSON data of LOD scan
#'
#' Each element in the array is:
#'
#' id: string (Ensembl ID)
#' chr: string
#' pos: integer
#' pheno1: numeric
#'
#' Example:
#'     {"data":[{"1_40055","1",40055,2.4787},
#'              {"1_87283","1",87283,2.4787},
#'              ...
#'              {"1_134512","1",134512,:2.4787}
#'             ]}
#'
#* @get /lodscan/protein
#* @post /lodscan/protein
http_lodscan_protein <- function(req, res, id, regress_local=FALSE, ncores=0) {
    # start the clock
    ptm <- proc.time()
    
    # get the Ensembl protein index 
    idx <- which(annot.protein$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }
    
    # make sure ncores is appropriate  
    num_cores = nvl_integer(ncores, 0)
    
    # set covars
    addcovar <- covar[,-1]
    
    # to regress local genotype, add neareast marker to covariates
    if (fix_boolean(regress_local)) {
        #addcovar <- cbind(addcovar, genoprobs[,-1,annot.protein$nearest_snp[idx]])
        mkr_tmp = as.character(snps[annot.protein$nearest_snp[idx],1])
        chr_tmp = as.character(snps[annot.protein$nearest_snp[idx],2])
        addcovar <- cbind(addcovar, probs[,chr_tmp][[1]][,-1,mkr_tmp])        
    }
    
    # perform the scan using QTL2
    temp <- (scan1(genoprobs=probs, 
                   kinship=Glist, 
                   pheno=expr.protein[,idx], 
                   addcovar=addcovar, 
                   cores=num_cores, 
                   reml=TRUE))
    
    # construct a 2 dimensional array of data
    temp_data <- data.table(id=snps$marker, chr=snps$chr, pos=snps$pos, temp)
    # setting colnames to NULL removes the names in the JSON and return an array
    colnames(temp_data) <- NULL
    to_return <- list(data=temp_data)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])
    
    return (to_return)
}

#' Get founder coefficients
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl gene identifier
#' @param chrom the chromosome
#' @param regress_local TRUE to regress local genotype
#' @param blup TRUE to perform Best Linear Unbiased Predictors 
#' @param center TRUE to center around 0
#' @param ncores number of cores to use (0=ALL)
#'
#' @return JSON data
#'
#' Each element in the array consists of:
#'
#' id: strieng (marker id)
#' chr: string
#' pos: numeric
#' A: numeric
#' B: numeric
#' C: numeric
#' D: numeric
#' E: numeric
#' F: numeric
#' G: numeric
#' H: numeric
#'
#' Example:
#' {"data":[{"id":"2_78746","chr":"2","pos":0.0787,"A":-0.0667,"B":0.1492,"C":0.1148,"D":0.1812,"E":-0.0165,"F":0.0199,"G":-0.2083,"H":-0.1736},
#'          ...
#'          {"id":"2_120469","chr":"2","pos":0.1205,"A":-0.0667,"B":0.1492,"C":0.1148,"D":0.1812,"E":-0.0165,"F":0.0199,"G":-0.2083,"H":-0.1736}
#'         ]}
#'
#* @get /foundercoefs
#* @post /foundercoefs
http_foundercoefs <- function(req, res, id, chrom, regress_local=FALSE, blup=FALSE, center=TRUE, ncores=0) {
    # start the clock
    ptm <- proc.time()

    # get the Ensembl gene index 
    idx <- which(annot.mrna$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }

    # make sure the chromosome data exists
    if (is.null(Glist[[chrom]])) {
        return (set_error(res, 400, paste0("chrom not found: ", chrom)))
    }

    # make sure ncores is appropriate  
    num_cores = nvl_integer(ncores, 0)
    
    # set covariates
    addcovar <- covar[,-1]
    
    # to regress local genotype, add neareast marker to covariates
    if (fix_boolean(regress_local)) {
        #addcovar <- cbind(addcovar, genoprobs[,-1,annot.mrna$nearest_snp[idx]])
        mkr_tmp = as.character(snps[annot.mrna$nearest_snp[idx],1])
        chr_tmp = as.character(snps[annot.mrna$nearest_snp[idx],2])
        addcovar <- cbind(addcovar, probs[,chr_tmp][,-1,mkr_tmp])
    }
    
    if (!fix_boolean(blup)) {
        temp <- scan1coef(genoprobs = probs[,chrom],
                          pheno = expr.mrna[,idx],
                          kinship=Glist[[chrom]],
                          addcovar = covar)
    } else {
        temp <- scan1blup(genoprobs = probs[,chrom],
                          pheno = expr.mrna[,idx],
                          kinship=Glist[[chrom]],
                          addcovar = covar,
                          cores=num_cores)
    }
    
    if (fix_boolean(center)) {
        temp[,LETTERS[1:8]] <- temp[,LETTERS[1:8]] - rowMeans(temp[,LETTERS[1:8]], na.rm=TRUE)
    }

    to_return <- list(data=data.table(id=names(map[[chrom]]), chr=chrom, pos=map[[chrom]], temp[,LETTERS[1:8]]))
    
    # stop the clock
    elapsed <- proc.time() - ptm
    print(paste0("Elapsed time: ", elapsed["elapsed"]))
    
    # return the data
    return (to_return)
}


#' Perform mediation analysis
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl gene identifier
#' @param mid a marker identifier
#'
#' return JSON data from mediation
#'
#' Example:
#'
#' [{"id":"ENSMUSG00000000001","symbol":"Gnai3","chr":"3","pos":108126713,"LOD":2.9505},
#'  ...
#'  {"id":"ENSMUSG00000000031","symbol":"H19","chr":"7","pos":142576836,"LOD":2.6687}
#' ]
#'
#* @get /mediate
#* @post /mediate
http_mediate <- function(req, res, id, mid) {
    # start the clock
    ptm <- proc.time()

    # get the Ensembl gene index 
    idx <- which(annot.mrna$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }

    # get the marker index 
    mrkx <- which(snps$marker == mid)
    chr_tmp = as.character(snps[mrkx,2])

    if (length(mrkx) == 0) {
        return (set_error(res, 400, paste0("mid not found: ", mid)))
    }

    annot <- annot.mrna[,c("id", "symbol", "chr")]
    annot$pos <- annot.mrna$middle_point

    
    # perform the mediation
    to_return <- mediation.scan(target=expr.mrna[,idx],  
                                mediator=expr.mrna, 
                                annotation=annot,
                                covar=covar, 
                                #qtl.geno=genoprobs[,,mrkx], 
                                qtl.geno=probs[[chr_tmp]][,,mid], 
                                verbose=FALSE)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])

    return (to_return)
}


#' Get the expression
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl gene identifier
#'
#' @return JSON data
#'
#' JSON data is dependent upon covariates.
#'
#' Example:
#'  {"data":[{"mouse_id":"F01","Sex":"F","Generation":"G4","Litter":2,"Diet":"hf","Coat.Color":"agouti","GenerationLitter":"G4_2","expression":0.2243,"_row":"F01"},
#'           ...
#'           {"mouse_id":"F02","Sex":"F","Generation":"G4","Litter":2,"Diet":"hf","Coat.Color":"black","GenerationLitter":"G4_2","expression":0.3498,"_row":"F02"}
#'          ]
#'
#* @get /expression/mrna
http_expression_mrna <- function(req, res, id) {
    # start the clock
    ptm <- proc.time()

    # get the index 
    idx <- which(annot.mrna$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }

    # simple to get the expression data and tack on
    output <- cbind(annot.samples, expression=expr.mrna[,idx])
    colnames(output)[1] <- ("mouse_id")

    # the types of expression data
    t <- list()
    for (f in covar_factors$column_name) {
        stopifnot(!is.null(annot.samples[[f]]))
        if (is.factor(annot.samples[[f]])) {
            t[[f]] <- levels(annot.samples[[f]])
        } else {
            t[[f]] <- unique(annot.samples[[f]])
        }
    }
    
    to_return <- list(data=output, data_types=t)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])
    
    return (to_return)
}


#' Get the expression
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl protein identifier
#'
#' @return JSON data
#'
#' JSON data is dependent upon covariates.
#'
#' Example:
#'  {"data":[{"mouse_id":"F01","Sex":"F","Generation":"G4","Litter":2,"Diet":"hf","Coat.Color":"agouti","GenerationLitter":"G4_2","expression":0.2243,"_row":"F01"},
#'           ...
#'           {"mouse_id":"F02","Sex":"F","Generation":"G4","Litter":2,"Diet":"hf","Coat.Color":"black","GenerationLitter":"G4_2","expression":0.3498,"_row":"F02"}
#'          ]
#'
#* @get /expression/protein
http_expression_protein <- function(req, res, id) {
    # start the clock
    ptm <- proc.time()
    
    # get the index 
    idx <- which(annot.protein$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }
    
    # simple to get the expression data and tack on
    output <- cbind(annot.samples, expression=expr.protein[,idx])
    colnames(output)[1] <- ("mouse_id")
    
    # the types of expression data
    t <- list()
    for (f in covar_factors$column_name) {
        stopifnot(!is.null(annot.samples[[f]]))
        if (is.factor(annot.samples[[f]])) {
            t[[f]] <- levels(annot.samples[[f]])
        } else {
            t[[f]] <- unique(annot.samples[[f]])
        }
    }
      
    to_return <- list(data=output, data_types=t)
    
    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])
    
    return (to_return)
}


#' Get the ids
#'
#' @param req the request object
#' @param res the response object
#' @param id_type 'mrna' or 'protein'
#'
#' @return JSON data
#'
#' Structure depends on the id_type.
#'
#' Example:
#' id_type = mrna 
#'     {"ids":[{"gene_id":"ENSMUSG00000000001"},{"gene_id":"ENSMUSG00000000031"},...{"gene_id":"ENSMUSG00000000049"}]}
#' id_type = protein 
#'     {"ids":[{"gene_id":"ENSMUSG00000000001", "protein_id":ENMUSP0000.."},...]}
#'
#* @get /ids
http_get_ids <- function(req, res, id_type) {
    # start the clock
    ptm <- proc.time()

    if (toupper(id_type) == "MRNA") {
        to_return <- (list(ids=list()))

        if (exists('annot.mrna')) {
            to_return <- (list(ids=data.frame(gene_id=annot.mrna$id)))
        }

        # stop the clock        
        elapsed <- proc.time() - ptm
        track_time(req, elapsed["elapsed"])
        
        return (to_return)
    } else if (toupper(id_type) == "PROTEIN") {
        to_return <- (list(ids=data.frame(protein_id=list(), gene_id=list())))

        if (exists('annot.protein')) {
            to_return <- (list(ids=data.frame(protein_id=annot.protein$id, gene_id=annot.protein$gene_id)))
        }
        
        # stop the clock
        elapsed <- proc.time() - ptm
        track_time(req, elapsed["elapsed"])

        return (to_return)
    } 
    
    return (set_error(res, 400, "Invalid id_type"))
}


#' Perform on SNP Association mapping (protein)
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl protein identifier
#' @param chrom the chromsome
#' @param location the location in base pairs
#' @param window_size how many base pairs (before and after) to perform scan
#' @param ncores number of cores to use (0=as many as there is)
#'
#' @return JSON data
#'
#' An array of objects where each object is:
#'     snp: string
#'     chr: string
#'     pos: decimal
#'     alleles: string
#'     sdp: integer
#'     ensembl_gene: string
#'     csq: string
#'     index: integer
#'     interval: integer
#'     on_map: boolean
#'     lod: decimal
#'
#' TODO: Trim the amount of data being passed back
#'
#' Example:
#' [{"snp":"rs33090065","chr":"5","pos":128.2978,"alleles":"A|G","sdp":233,"ensembl_gene":"ENSMUSG00000034310","csq":"intron_variant","index":1,"interval":2823,"on_map":true,"lod":0.3721},
#'  ...
#'  {"snp":"rs266155287","chr":"5","pos":128.2978,"alleles":"C|T","sdp":128,"ensembl_gene":"ENSMUSG00000034310","csq":"intron_variant","index":2,"interval":2823,"on_map":false,"lod":0.3026}
#' ]
#'
#* @get /snpassoc/protein
#* @post /snpassoc/protein
http_snp_assoc_mappin_protein <- function(req, res, id, chrom, location, window_size=500000, ncores=0) {

    # start the clock
    ptm <- proc.time()

    # protein id
    idx <- which(annot.protein$id == id)

    # replace with gene id
    idx <- which(annot.mrna$id == annot.protein[idx,]$gene_id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }
    
    sel.chr <- chrom
    
    loc <- nvl_integer(location, -1)
    if (loc == -1) {
        return (set_error(res, 400, paste0("location is invalid: ", location)))
    }
    pos <- loc / 1000000.0

    windowsize <- nvl_integer(window_size, -1)
    if (windowsize == -1) {
        return (set_error(res, 400, paste0("window_size is invalid: ", window_size)))
    }

    window.length <- windowsize/1000000.0
    window.range <- pos + c(-1,1)*window.length

    # make sure ncores is appropriate  
    num_cores = nvl_integer(ncores, 0)

    # extract SNPs from the database
    my_db <- src_sqlite(db.file, create = FALSE)
    window.snps <- tbl(my_db, sql("SELECT * FROM snps")) %>%
        filter(chr==sel.chr, pos_Mbp>=window.range[1], pos_Mbp<=window.range[2]) %>%
        arrange(pos_Mbp) %>%
        collect(n=Inf)

    colnames(window.snps)[c(1,3)] = c("snp", "pos")
    window.snps = index_snps(map = map, window.snps)
    #window.snps = as.data.frame(window.snps)

    addcovar <- covar[,-1]

    # to regress local genotype, add neareast marker to covariates
    # if (regress_local) addcovar <- cbind(addcovar, genoprobs[,-1,annot.mrna$middle_point[idx]])

    # convert allele probs to SNP probs
    snppr <- genoprob_to_snpprob(probs, window.snps)

    # finally the scan
    out_snps <- scan1(pheno=expr.mrna[,idx], kinship = Glist[[sel.chr]], genoprobs=snppr, addcovar=addcovar, cores=ncores)

    map_tmp <- qtl2scan:::snpinfo_to_map(window.snps)
    tmp <- qtl2plot:::expand_snp_results(out_snps, map_tmp, window.snps)

    to_return <- window.snps
    to_return$lod <- tmp$lod[,1]

    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])
    
    return (to_return)
    
}

#' Perform on SNP Association mapping (protein)
#'
#' @param req the request object
#' @param res the response object
#' @param id an Ensembl protein identifier
#' @param chrom the chromsome
#' @param location the location in base pairs
#' @param window_size how many base pairs (before and after) to perform scan
#' @param ncores number of cores to use (0=as many as there is)
#'
#' @return JSON data
#'
#' An array of objects where each object is:
#'     snp: string
#'     chr: string
#'     pos: decimal
#'     alleles: string
#'     sdp: integer
#'     ensembl_gene: string
#'     csq: string
#'     index: integer
#'     interval: integer
#'     on_map: boolean
#'     lod: decimal
#'
#' TODO: Trim the amount of data being passed back
#'
#' Example:
#' [{"snp":"rs33090065","chr":"5","pos":128.2978,"alleles":"A|G","sdp":233,"ensembl_gene":"ENSMUSG00000034310","csq":"intron_variant","index":1,"interval":2823,"on_map":true,"lod":0.3721},
#'  ...
#'  {"snp":"rs266155287","chr":"5","pos":128.2978,"alleles":"C|T","sdp":128,"ensembl_gene":"ENSMUSG00000034310","csq":"intron_variant","index":2,"interval":2823,"on_map":false,"lod":0.3026}
#' ]
#'
#* @get /snpassoc/mrna
#* @post /snpassoc/mrna
http_snp_assoc_mapping_mrna <- function(req, res, id, chrom, location, window_size=500000, ncores=0) {

    # start the clock
    ptm <- proc.time()

    idx <- which(annot.mrna$id == id)
    
    if (length(idx) == 0) {
        return (set_error(res, 400, paste0("id not found: ", id)))
    }
    
    sel.chr <- chrom
    
    loc <- nvl_integer(location, -1)
    if (loc == -1) {
        return (set_error(res, 400, paste0("location is invalid: ", location)))
    }
    pos <- loc / 1000000.0

    windowsize <- nvl_integer(window_size, -1)
    if (windowsize == -1) {
        return (set_error(res, 400, paste0("window_size is invalid: ", window_size)))
    }

    window.length <- windowsize/1000000.0
    window.range <- pos + c(-1,1)*window.length

    # make sure ncores is appropriate  
    num_cores = nvl_integer(ncores, 0)

    # extract SNPs from the database
    my_db <- src_sqlite(db.file, create = FALSE)
    window.snps <- tbl(my_db, sql("SELECT * FROM snps")) %>%
        filter(chr==sel.chr, pos_Mbp>=window.range[1], pos_Mbp<=window.range[2]) %>%
        arrange(pos_Mbp) %>%
        collect(n=Inf)

    colnames(window.snps)[c(1,3)] = c("snp", "pos")
    window.snps = index_snps(map = map, window.snps)
    #window.snps = as.data.frame(window.snps)

    addcovar <- covar[,-1]

    # to regress local genotype, add neareast marker to covariates
    # if (regress_local) addcovar <- cbind(addcovar, genoprobs[,-1,annot.mrna$middle_point[idx]])

    # convert allele probs to SNP probs
    snppr <- genoprob_to_snpprob(probs, window.snps)

    # finally the scan
    out_snps <- scan1(pheno=expr.mrna[,idx], kinship = Glist[[sel.chr]], genoprobs=snppr, addcovar=addcovar, cores=ncores)

    map_tmp <- qtl2scan:::snpinfo_to_map(window.snps)
    tmp <- qtl2plot:::expand_snp_results(out_snps, map_tmp, window.snps)

    to_return <- window.snps
    to_return$lod <- tmp$lod[,1]

    # stop the clock
    elapsed <- proc.time() - ptm
    track_time(req, elapsed["elapsed"])
    
    return (to_return)
    
}

#' Get the phenotypes
#'
#' @param req the request object
#' @param res the response object
#'
#' NOT COMPLETED
#'
#' @return phenotypes
#'
#* @get /phenotypes
http_get_phenotypes <- function(req, res) {
    # start the clock
    ptm <- proc.time()

    if (exists('annot.pheno')) {
        to_return <- annot.pheno[which(annot.pheno$omit == FALSE),]

        # stop the clock        
        elapsed <- proc.time() - ptm
        track_time(req, elapsed["elapsed"])
        
        return (to_return)
    }

    return (set_error(res, 400, "Invalid id_type"))
}


#
# any duplicated gene ids in protein data
#
# any(duplicated(annot.protein$gene_id))
#

# list all duplicated protein gene ids
#
# annot.protein$gene_id[which(duplicated(annot.protein$gene_id)==TRUE)]
#

#
# list gene ids that are IN protein data but NOT IN mrna data
#
# annot.protein$gene_id[which((annot.protein$gene_id %in% annot.mrna$id)==FALSE)]
#

#
# list gene ids that are IN mrna data but NOT IN protein data
#
# annot.mrna$id[which((annot.mrna$id %in% annot.protein$gene_id)==FALSE)]
#
