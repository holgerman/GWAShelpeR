#' @title CALCULATE ASSOCIATION WITH PLINK 2.0 FROM PGEN FILE
#' @description Uses a psam file and plink 2.0 to quickly calculate associations for quantitative and binary outcome variable controlling for covariates for autosomes and the X-chromosome
#' @param sampleinfo_fn  File with outcome variable and covariates mit phenotypen, first two columns used as ID, must include the same IDs as defined in the psam file specified in `pfiles` parameter
#' @param covariables colnames of sampleinfo_fn specifying confounder to be included in model
#' @param adjust4sex Shall be adjusted for Sex? Do not include Sex as covariables, as this might result in bad behaviour in X-chromosomal analysis, Default: T
#' @param outcome colnames of sampleinfo_fn specifying the outcome variable, can be binary or continous
#' @param out_fn directory and prefix of result filenames
#' @param keep_fn Filename of a file consisting of two columns used as ID of individuals to keep. `''` means no filtering, Default: ''
#' @param snps2extract Filename of a file including SNP IDs as specified in the pvar file specified in parameter `pfiles` to keep. `''` means no filtering, Default: ''
#' @param pfiles full path to psam, pgen, and pvar fileset including the filename without file extension, Default: "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1612_lifea1_genotypisierungsrunde_3/Imputation/03_imputed_prephased/bgen1_2/combined/s302_3_PLINK_QCed_StandardRS_combined_a1_rd3_chr1to23_N7660"
#' @param unixcomputer name of compute server, required to choose right plink version. Can be 'amanMRO' that will lead to use of version specified in  plinklocation_avx, otherwise, plink version specified in plinklocation_unix64 is used, Default: 'amanMRO'
#' @param use_firth_fallback to stabilize low maf calculation in logistic regression, see https://www.cog-genomics.org/plink/2.0/assoc, Default: T
#' @param code_chrX encoding can be either `m:0_2_f:0_1_2` or `m:0_1_f:0_1_2`, Default: 'm:0_2_f:0_1_2'
#' @param showCovarStats Include association summary statistic of covariates in output file, Default: F
#' @param max_vif variance inflation factor still ok in plink regression, reflecting correlation of covariates, can be calculated for a model using car::vif(mymodel) , Default: 50
#' @param morePlinkParameter plink parameters added to the call, see https://www.cog-genomics.org/plink/2.0/, Default: ' --threads 2 --memory 6000'
#' @param gz compress result files to gz, Default: T
#' @param plinklocation_avx Location of plink executable requiring unix 64 bit architecture with avx support, Default: '/net/ifs1/san_projekte/projekte/genstat/07_programme/plink2.0/20201105/unix_64/plink2'
#' @param plinklocation_unix64 Location of plink executable requiring any unix 64 bit architecture, Default: '/net/ifs1/san_projekte/projekte/genstat/07_programme/plink2.0/20201105/unix_avx2/plink2'
#' @return Returns a data frame with standard plink association results, see also https://www.cog-genomics.org/plink/2.0/formats#glm_linear or https://www.cog-genomics.org/plink/2.0/formats#glm_logistic Additionally all plink output is saved in the location specified in `out_fn`
#' @details Does not use rounded genotypes, but imputed gene doses. Works currently only on unix computeserver of GEnstat. Function in progress. Last modificaiton: 2020-11-11
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calcPlink2AssocFromPGEN
#' @export
#' @importFrom testthat test_that expect_equal
#' @importFrom data.table fread fwrite
#' @importFrom toolboxH formateTimediff moveColFront
#' @importFrom stringr str_split str_locate_all str_sub str_replace_all

#
calcPlink2AssocFromPGEN <- function(sampleinfo_fn,
                                    covariables ,
                                    adjust4sex  = T,
                                    outcome ,
                                    out_fn,
                                    keep_fn="",
                                    snps2extract = "",
                                    pfiles = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1612_lifea1_genotypisierungsrunde_3/Imputation/03_imputed_prephased/bgen1_2/combined/s302_3_PLINK_QCed_StandardRS_combined_a1_rd3_chr1to23_N7660",
                                    unixcomputer = "amanMRO",
                                    use_firth_fallback = T,
                                    code_chrX = "m:0_2_f:0_1_2",
                                    showCovarStats = F,
                                    max_vif = 50 , # default 50
                                    morePlinkParameter = " --threads 2 --memory 6000",
                                    gz = T,
                                    plinklocation_avx = "/net/ifs1/san_projekte/projekte/genstat/07_programme/plink2.0/20201105/unix_64/plink2",
                                    plinklocation_unix64 = "/net/ifs1/san_projekte/projekte/genstat/07_programme/plink2.0/20201105/unix_avx2/plink2"
) {

  time1 = Sys.time()
  psamfile_fn = paste0(pfiles, ".psam")

  if(unixcomputer == "amanMRO") callplink = plinklocation_avx  else callplink = plinklocation_unix64
  message("Using plink ", callplink)


  check = system(callplink, show.output.on.console = F)
  testthat::test_that("Plink program starts up", testthat::expect_equal(check, 1))


  if(use_firth_fallback==T) use_firth_fallback_text ="firth-fallback"  else use_firth_fallback_text =  "no-firth"


  if(code_chrX %in% c( "m:0_2_f:0_1_2")) code_chrX_text = "--xchr-model 2"  else  if(code_chrX == "m:0_1_f:0_1_2") code_chrX_text = "--xchr-model 1" else stop('Analysismode plink accepts only values  "m:0_2_f:0_1_2", and "m:0_1_f:0_1_2" for paramter "code_chrX"')

  if(length(snps2extract) ==0) extract_text = "" else {
    goodsnps_fn = paste0(out_fn, ".snps2extract")
    extract_text = paste("--extract", goodsnps_fn)
    fwrite(data.table(snp=snps2extract), goodsnps_fn, col.names = F)
  }

  if(keep_fn =="") keeptext = "" else keeptext = paste("--keep", keep_fn)


  sampleinfo = data.table::fread(sampleinfo_fn)
  setnames(sampleinfo, names(sampleinfo)[1:2], c('FID', 'IID'))

  psamfile  = data.table::fread(psamfile_fn)
  setnames(psamfile, names(psamfile)[1:2], c('FID', 'IID'))

  # qlist1 = venn2(sampleinfo[,paste(get(names(sampleinfo)[1]),get(names(sampleinfo)[2]))], psamfile[,paste(get(names(psamfile)[1]),get(names(psamfile)[2]))], mylabels = c("sampleinfo IDs", 'psamfile IDs'), plotte=F)
  #
  # testthat::test_that("Same prim.&sec.IDs in sampleinfo and psam file", testthat::expect_equal(length(c(qlist1$q2, qlist1$q3)), 0))

  setDF(psamfile)
  setDF(sampleinfo)

  pheno = sampleinfo[,c('FID', 'IID', outcome)]

  if(keep_fn != "") {
    keep = data.table::fread(keep_fn, header = F)
    setDF(keep)
    keep$id = paste(keep[,1], keep[,2])
    psamfile = psamfile[paste(psamfile[,"FID"], psamfile[,"IID"]) %in% keep$id,]
    pheno = pheno[paste(pheno[,"FID"], pheno[,"IID"]) %in% keep$id,]
    # if(covar_fn !="") covar = covar[paste(covar[,1], covar[,2]) %in% keep$id,]
    #
  }


  qlist1= venn2(paste(psamfile[,"FID"], psamfile[,"IID"]), paste(pheno[,"FID"], pheno[,"IID"]) , plotte = F)

  if(length(c(qlist1$q2, qlist1$q3)) >0) {
    print(str(qlist1))
    stop("check that the same individuals with the same pedigree IDs are in files\n",paste(c(sample_fn, sampleinfo_fn), collapse = "\n"))
  }


  if(identical(covariables ,"")) covartext ="" else {
    stopifnot(all(covariables %in% names(sampleinfo)))
    if(adjust4sex==F) message("using covars ", paste(covariables, collapse = " + "))
    if(adjust4sex==T) message("using covars ", paste(c("[sex from samplefile]",covariables), collapse = " + "))

    covar_fn = paste0(out_fn,".covars")
    covartext = paste("--covar", covar_fn)
    covar = sampleinfo[, c("FID", "IID", covariables)]

    write.delim(covar, covar_fn)
  }

  message("using outcomes\n", paste(outcome, collapse = "\n"))

  pheno_fn = paste0(out_fn, ".phenos")
  phenotext = paste("--pheno", pheno_fn)
  pheno = sampleinfo[, c(1,2, which(names(sampleinfo) %in% outcome))]
  write.delim(pheno, pheno_fn)



  sextext = if(adjust4sex==T) sextest = "sex" else sextext = ''



  if(showCovarStats) showCovarStats_text = "" else showCovarStats_text = "hide-covar"

  #call2----
  call2 = paste(callplink, '--pfile',psamfile_fn = pfiles, "--glm ",sextext," ", showCovarStats_text,"cols=chrom,pos,ref,alt,firth,test,nobs,machr2,a1freq,a1freqcc,a1countcc,orbeta,se,ci,tz,p",use_firth_fallback_text,code_chrX_text, " --ci 0.95 --pheno",pheno_fn,covartext,keeptext,extract_text, morePlinkParameter,"--vif ",max_vif," --covar-variance-standardize --out", out_fn)


  message("-----------\nRunnning Plink command:\n",call2)


  system(call2)


  message("---------------------------------\nTotal time for calculating genedose-associations :\n", toolboxH::formateTimediff(Sys.time()-time1))

  outpattern = stringr::str_split(out_fn,"/") %>% unlist %>% last

  allresi = lapply(outcome, function(myoutcome) {
    # myoutcome = outcome[1]
    message("collecting results from ", myoutcome)
    myoutpattern = paste0(outpattern, ".", myoutcome)
    message("Using grep-pattern to identify plink resultfile:\n",myoutpattern) # myoutpattern = 's302_1_candidatelookup'

    position_slash = stringr::str_locate_all(out_fn, "/")
    lastposition_slash = max(position_slash[[1]][,"end"])
    out_fn_tocheck = stringr::str_sub(out_fn, 1, lastposition_slash)
    my_resultfile_pre = dir(path = out_fn_tocheck, recursive = T, pattern = myoutpattern)
    my_resultfile = paste0(out_fn_tocheck, "/", my_resultfile_pre)

    message("Checking guessed name of  PLINK outfile,:\n",my_resultfile)

    if(length(my_resultfile)==2) {
      my_resultfile = stringr::str_replace_all(my_resultfile, "\\.gz", "") %>% unique

    }
    if(file.exists(my_resultfile)==F) {


      stop("Did not identified correctly name of PLINK outfile, as the guessed filename does not exists:\n",my_resultfile)
    }

    resi = data.table::fread(my_resultfile)

    if(gz==T & (grepl("$\\.gz", my_resultfile)==F)) {
      system(paste("gzip -f", my_resultfile))
      my_resultfile = paste0(my_resultfile, ".gz")
    }
    resi$out_fn = my_resultfile
    resi$pheno = myoutcome
    resi = toolboxH::moveColFront(resi, "pheno")
    resi
  }) %>% rbindlist(., use.names = T)


  allresi

  ## checke ob alle guten snps da sind
  snps2extract_dt = data.table(snps = snps2extract)
  allresi2 = merge(allresi, snps2extract_dt, by.x  = "ID", by.y = "snps", all = T, sort = F)


  out_fn_results = paste0(out_fn, '_ALLPHENOS.txt')
  message("\n---------------------\nSaving all results in one file here:\n", out_fn_results)
  data.table::fwrite(allresi2, out_fn_results, sep = "\t")

  allresi2

  if(gz ==T)      {
    message("GZ ing all result file...")
    system(paste("gzip -f", out_fn_results))

  }

  allresi2
}

