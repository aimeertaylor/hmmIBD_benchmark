plot_both <- function(to_file) {
    # some hardwired sample numbers
#    npair_cam = 79401
#    npair_ghana = 31375
    npair_cam_ghana = 100945

    chr = 7
    if (to_file == 1) {
      outfile = "cam_ghana_1v2.pdf"
      pdf(outfile);
    }
    gene = c(403222,406317) / 1000
    xlim = c(0,750)
    tag = "pfcrt"
    height = c(0,0)


    file = "twopop.overlaps"
    d = read.delim(file)
    starts2 = d$local_start[d$chrom == chr]/1000
    ends2 = d$local_stop[d$chrom == chr]/1000
    n_ident2 = d$N_ident_pair[d$chrom == chr] / npair_cam_ghana
    
    file = "onepop.overlaps"
    d = read.delim(file)
    starts1 = d$local_start[d$chrom == chr]/1000
    ends1 = d$local_stop[d$chrom == chr]/1000
    n_ident1 = d$N_ident_pair[d$chrom == chr]  / npair_cam_ghana

    ymax = max(n_ident2, n_ident1, na.rm=TRUE)
    ylim = c(0, ymax)

    nseg = length(starts2)
    nidentm=rbind(n_ident2,n_ident2)
    posm = rbind(starts2, ends2)
    pos_line = as.vector(posm)
    nident_line = as.vector(nidentm)
    matplot(posm,nidentm,col="blue",type="l",ylim=ylim, lty=1,
    ylab="Fraction of sample pairs IBD", xlab=paste("Position on chromosome",chr,"(kb)"), 
    main=paste("Fraction of pairs that are IBD,",tag,"region"), 
    lwd=2,
    xlim=xlim)
    points(pos_line, nident_line, type="l", col="black",lwd=.5)

    nidentm = rbind(n_ident1, n_ident1)
    posH = rbind(starts1, ends1)
    matplot(posH,nidentm,col="red",type="l",add=TRUE, lty=1,
    xaxt="n",lwd=2)
    pos_line = as.vector(posH)
    nident_line = as.vector(nidentm)
    points(pos_line, nident_line, type="l", col="red",lwd=.5)

    points(gene,height,col="green",type="l",lwd=4)

    legend("topright", c("two populations", "one population", tag),
        col=c("blue","red", "green"), lty=c(1,1,1), lwd=c(2))

    if(to_file == 1) {dev.off();}
}
