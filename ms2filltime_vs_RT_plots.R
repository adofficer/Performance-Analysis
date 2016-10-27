require(ggplot2)
require(data.table)
require(doBy)
require(RColorBrewer)
require(hexbin)

"
Small script that plots the distributions of MS2 fill time vs RT
Includes the median values for each RT subset as well

Assumptions: input file contains numeric columns that are named
             tripFillMsec and retetionTimeMin.

Input: Ideally a Spec Features file output from Spectrum Mill, several
       optional arguments to specify display options.  Can take any file
       containing numeric columns with names trapFillMsec and retentionTimeMin

Output: displays a user-defined number of histograms (changed by cutsize)
        that correspond to the MS2 fill time distributions within each
        RT bin
"


plotMS2filltimevsRTbins <- function(specfile = "SpecFeatures.1.tsv", hitfile = "hitTable.1.tsv",
                                    cutsize = 10, gradlength = 110) {
  
  # read the data in
  hitspecfeat <- readSpecandHitFiles(specfile, hitfile)
  
  # sets up cuts along the RT dimension (retbins)
  retcuts <- seq(0, gradlength, by = cutsize)
  # dealing with the situation where gradlength %% cutsize != 0
  if (max(retcuts) != gradlength) {
    retcuts <- c(retcuts, gradlength)
  }
  cutvector <- cut(hitspecfeat$retentionTimeMin,
                  breaks = retcuts)
  hitspecfeat$retbins <- cutvector

  # computes median MS2 fill time in each retention time bin
  retbinmedians <- summaryBy(trapFillMsec ~ retbins,
                            data = hitspecfeat,
                            FUN = median)

  # plot results
  ms2filltimehist <- ggplot(hitspecfeat,
                            aes(x = trapFillMsec, fill = retbins))

  ms2filltimehist + geom_histogram(binwidth = (gradlength / 25)) +
                  facet_wrap(~retbins) +
                  geom_vline(data = retbinmedians, 
                              aes(xintercept = retbinmedians$trapFillMsec.median))
}

hexbinMS2filltimevsRT <- function(specfile = "SpecFeatures.1.tsv", hitfile = "hitTable.1.tsv",
                                  gradlength = 110, usehitfile = T) {
  # initialize the color spectrum for plot density
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  colorscheme <- rf(32)

  # read the file in
  hitspecfeat <- readSpecandHitFiles(specfile, hitfile, usehitfile)
  
  # bin the data
  hex_hits <- hexbin(hitspecfeat$retentionTimeMin, hitspecfeat$trapFillMsec, xbnds = c(0, gradlength),
                     ybnds = c(0, max(hitspecfeat$trapFillMsec)))

  # visualize the data
  plot(hex_hits, colramp = rf, colorcut = 15)
}

"
#setup const bin scales - round up to nearest 100
powerOfTenBase<-floor(log(max(h@count))/log(10))
prelimRounded<-signif(max(h@count),digits = (powerOfTenBase-1))
if (prelimRounded < max(h@count)) { prelimRounded<-prelimRounded+100}

H<-data.frame(cell=h@cell,prec_count=h@count)
J<-data.frame(cell=j@cell,hit_count=j@count)
K<-merge(H,J,all.x=TRUE)
K$hit_count[is.na(K$hit_count)]<-0
K$pct<-round(K$hit_count/K$prec_count*100)
k@count<-K$pct

plot(h,colramp=rf,xlab='RT (min)',ylab='precursor m/z', main=paste('Precursors Sampled',getwd(),sep=),mincnt=0,maxcnt=prelimRounded,colorcut=11)
"

readSpecandHitFiles <- function(specfile, hitfile, usehitfile = T) {
  
  # checks if the file exists, hopefully only error handling I have to do :3
  if (file.exists(specfile)) {
    specfeat <- fread(specfile, 
                      select = c("trapFillMsec", "retentionTimeMin", "originalSpectrumFile"))
  } else {
    stop("The spec file does not exist!")
  }
  
  if (usehitfile) {
    if (file.exists(hitfile)) {
      hitlist <- fread(hitfile)
    } else {
      stop("The hit file does not exist!")
    }
    
  }
  
  if (usehitfile) {
    # string munging to match for the subsetting action
    specfeat$originalSpectrumFile <- paste(specfeat$originalSpectrumFile, ".spo", sep = "")
    
    # filter out only those spectra that were valid hits
    specfeat <- subset(specfeat, originalSpectrumFile %in% hitlist$mstagOutputFile)
  }
  
  return(specfeat)
}
