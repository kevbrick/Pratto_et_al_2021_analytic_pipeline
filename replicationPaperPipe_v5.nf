#!/usr/bin/env nextflow

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "PRATTO, BRICK ET AL. 2020 PIPELINE (Version 5.0)                          "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run replicationPaperPipe_v5.nf "
  log.info "             -c accessoryFiles/config/nextflow.config.nf "
  log.info "             -profile singularity "
  log.info "             --outdir output "
  log.info "             --projectdir `pwd` "
  log.info "             -with-trace -with-timeline -with-report"
  log.info " "
  log.info "HELP: nextflow run replicationPaperPipe_v5.nf --help"
  log.info " "
  log.info "================================================================================================================="
  log.info "Required Arguments:"
  log.info " --projectdir    Top-level folder (contains /accessoryFiles folder; default = .)"
  log.info " -profile modules / singularity"
  log.info " "
  log.info "Output Arguments"
  log.info "  --outdir       Output folder "
  log.info " "
  log.info "================================================================================================================"
  exit 1
  }

// Define arg defaults
params.threads        = 6
params.mem            = "16G"
params.debugmode      = true

//output and tempory directories
params.projectdir        = "."
params.tmpdir            = "/lscratch/\$SLURM_JOBID"

params.accessorydir      = "${params.projectdir}/accessoryFiles/"
params.mmOriSSDSdir      = "${params.accessorydir}/data/oriSSDS/mouse/"
params.RTdir             = "${params.accessorydir}/RTSeq/"
params.RDorgdir          = "${params.accessorydir}/RTSeq/replicationdomain.org"

params.codedir           = "accessoryFiles/scripts/"
params.codeRdir          = "accessoryFiles/scripts/R"
params.datadir           = "accessoryFiles/data/"
params.genomedir         = "accessoryFiles/genomeFiles/"
params.annotationdir     = "accessoryFiles/annotation/"
params.rnadir            = "accessoryFiles/RNASeq/"
params.modelsdir         = "accessoryFiles/RTSeq/bestModels/"
params.dataRTdir         = "accessoryFiles/RTSeq/"

Channel.fromPath("${params.accessorydir}")
       .into {af1; af2 ; af3 ; af4 ; af5 ; af6 ; af7 ; af8 ; af9 ; af10;
             af11; af12; af13; af14; af15; af16; af17; af18; af19; af20;
             af21; af22; af23; af24; af25; af26; af27; af28; af29; af30;
             af31; af32; af33; af34; af35; af36; af37; af38; af39; af40;
             af41; af42; af43; af44; af45; af46; af47; af48; af49; af50}

params.outdir            = ""
params.outdirPeaks       = "${params.outdir}/peakCalls"
params.outdirOrigins     = "${params.outdir}/origins"
params.outdirAnnot       = "${params.outdir}/annotation"
params.outdirHiC         = "${params.outdir}/annotation/hic"
params.outdirGomez       = "${params.outdirAnnot}/gomezOrigins"
params.outdirBW          = "${params.outdir}/bigwig"
params.outdirBG          = "${params.outdir}/bedgraph"
params.outdirRDBG        = "${params.outdir}/bedgraph/replicationDomain.org"
params.outdirModelBG     = "${params.outdir}/bedgraph/RTsim"
params.outdirSlices      = "${params.outdir}/slices"
params.outdirDTdata      = "${params.outdir}/deeptoolsData"
params.outdirRTables     = "${params.outdir}/Rtables"
params.outdirModel       = "${params.outdir}/model"
params.outdirWGGIF       = "${params.outdir}/model/gif/WG"
params.outdirCSGIF       = "${params.outdir}/model/gif/CS"
params.outdirImages      = "${params.outdir}/figs"
params.outdirFigs        = "${params.outdir}/figs"
params.outdirModelFigs   = "${params.outdir}/figs/RTsim"
params.outbam            = "${params.outdir}/bam"

params.reps              = 1
params.test              = false
params.genome            = 'mm10'

params.outdir_tmp        = "/tmp"

params.getRNASeqData     = true
params.getStatsAgain     = false
params.getHumanDMC1      = false
params.callHumanHS       = false

def mm10FA               = "${params.genomedir}/mm10_genome.fa"
def mm10IDX              = "${params.genomedir}/mm10_genome.fa.fai"
def hg38FA               = "${params.genomedir}/hg38_genome.fa"
def hg38IDX              = "${params.genomedir}/hg38_genome.fa.fai"

//output and tempory directories
log.info "===================================================================="
log.info "Pratto, Brick et al. DNA replication in meiosis pipeline : "
log.info "===================================================================="
log.info "Sample exectution : "
log.info " "
log.info "nextflow run replicationPaperPipe_v4.nf "
log.info "             -c accessoryFiles/config/nextflow.config.nf"
log.info "             -profile singularity "
log.info "             --projectdir `pwd` "
log.info " "
log.info "===================================================================="
log.info "- nextflow args ----------------------------------------------------"
log.info "Nextflow config file : ${workflow.configFiles}"
log.info " "
log.info "- pipeline args ----------------------------------------------------"
log.info "Project dir          : ${params.projectdir}"
log.info "threads              : ${params.threads}"
log.info "mem                  : ${params.mem}"
log.info " "
log.info "- annotation data --------------------------------------------------"
log.info "code dir             : ${params.codedir}"
log.info "data dir             : ${params.datadir}"
log.info "genome dir           : ${params.genomedir}"
log.info "RNA-Seq dir          : ${params.rnadir}"
log.info " "
log.info "- output directories -----------------------------------------------"
log.info "outdir               : ${params.outdir}"
log.info "outdir [Peaks]       : ${params.outdirPeaks}"
log.info "outdir [Annotation]  : ${params.outdirAnnot}"
log.info "outdir [BG and BW]   : ${params.outdirBW}"
log.info "outdir [Images]      : ${params.outdirImages}"
log.info " "
log.info "--------------------------------------------------------------------"

if (params.getStatsAgain){
  log.info " "
  log.info "*********************************************************************"
  log.info "***************  WARNING !!! WARNING   ******************************"
  log.info " "
  log.info "getStatsAgain Enabled !!! This will cause a long runtime (+50hr) ....."
  log.info " "
  log.info "*********************************************************************"
  log.info "*********************************************************************"
  log.info " "
  }

// Build BED / BAM input channels
// Channel of oriSSDS BAMs
Channel
  .fromPath("${params.mmOriSSDSdir}/*.ssDNA_type1.bam")
  .ifEmpty { exit 1, "mm10 oriSSDS BAM files not found in ${params.mmOriSSDSdir} or mis-named.\nTry running getOriSSDSData.sh in the accessoryFiles/data/oriSSDS folder" }
  .map { sample -> tuple(getName(sample), 'mm10', file(getIDX(sample)), sample) }
  .into {oriSSDSBAMs_mm10_1; oriSSDSBAMs_mm10_2; ssdsBAMs; ssdsBAMs_b}

  def getName( file ) {
    def nm="${file.name}"
    def nRet = nm.replaceFirst(".oriSSDS.ssDNA_type1.bam","")
    return nRet
  }

  def getIDX( file ) {
    def nm="${file}"
    return nm.replaceFirst(/.bam/,".bam.bai")
  }

// Channel of oriSSDS BEDs
Channel
  .fromPath("${params.mmOriSSDSdir}/*.ssDNA_type1.bed")
  .ifEmpty { exit 1, "mm10 oriSSDS BED files not found or mis-named. Try running getOriSSDSData.sh in the accessoryFiles/data/oriSSDS folder" }
  .map { sample -> tuple(getBEDName(sample), 'mm10', sample) }
  .into {oriSSDSBEDs_mm10_1; oriSSDSBEDs_mm10_2; ssdsBEDs; ssdsBEDs_a}

  def getBEDName( file ) {
    def nm="${file.name}"
    def nRet = nm.replaceFirst(".oriSSDS.ssDNA_type1.bed","")
    return nRet
  }

// Channel of oriSSDS BEDs (used for slices)
Channel
  .fromPath("${params.mmOriSSDSdir}/*.ssDNA_type1.bed")
  .ifEmpty { exit 1, "mm10 oriSSDS BED files not found or mis-named. Try running getOriSSDSData.sh in the accessoryFiles/data/oriSSDS folder" }
  .set {oriSSDSBEDs_mm10_4Slices}

// Channel of oriSSDS BED folders
Channel
  .from( params.mmOriSSDSdir )
  .ifEmpty { exit 1, "oriSSDS BED folders not found or mis-named. Try running getOriSSDSData.sh in the accessoryFiles/data/oriSSDS folder" }
  .map { sample -> tuple(sample, 'mm10', 'mm10_OriSSDS') }
  .set {oriSSDSPaths}

Channel
  .fromPath("${params.accessorydir}/RTSeq/bedgraph/final/forModel/*bedgraph")
  .ifEmpty { exit 1, "RT bedgraphs NOT found" }
  .set {rtBGs}

Channel
  .fromPath("${params.RDorgdir}/mm10/*bedgraph")
  .ifEmpty { exit 1, "RD.org bedgraphs NOT found [${params.RDorgdir}/mm10/*bedgraph]" }
  .set {rdOrgBG_mm10}

Channel
  .fromPath("${params.RDorgdir}/hg38/*bedgraph")
  .ifEmpty { exit 1, "RD.org bedgraphs NOT found [${params.RDorgdir}/hg38/*bedgraph]" }
  .set {rdOrgBG_hg38}

// Set saturation curve params for origin calling
if (params.test =~ /test/){
  satCurvePCs = Channel.from(0.01,1.00)
  params.reps = 1
}else{
  satCurvePCs = Channel.from(0.01,0.025,0.05,0.075,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00)
  satCurvePCs = Channel.from(0.01,0.05,1.00)
}

process makeWinFiles {
  input:
  file (af) from af1.collect()

  output:
  file "win*bed" into (winFiles_a, winFiles_b, winFiles_c, winFiles_d, winFiles_e, winFiles_f)

  script:
  """
  #bedtools makewindows -g ${mm10IDX} -w 25          | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 25  )\' >win25.bed
  #bedtools makewindows -g ${mm10IDX} -w 200         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 200 )\' >win200.bed

  bedtools makewindows -g ${mm10IDX} -w 147         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 147 )\' >win147.bed
  bedtools makewindows -g ${mm10IDX} -w 500         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500 )\' >win500.bed
  bedtools makewindows -g ${mm10IDX} -w 500 -s 50   | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500 )\' >win500s50.bed
  bedtools makewindows -g ${mm10IDX} -w 1000 -s 100 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000)\' >win1ks100.bed
  bedtools makewindows -g ${mm10IDX} -w 1000 -s 147 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000)\' >win1ks147.bed
  bedtools makewindows -g ${mm10IDX} -w 2000 -s 100 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 2000)\' >win2ks100.bed
  """
  }

process getOtherOrigins {

  publishDir params.outdirGomez, mode: 'copy', overwrite: true

  input:

  output:
  file 'ori_ESC_gomez.bedgraph'    into gomezOriES
  file 'ori_MEF_gomez.bedgraph'    into gomezOriMEF
  file 'AlmeidaESC.peaks.mm10.bed' into ori_AlmeidaESC
  file 'AlmeidaMEF.peaks.mm10.bed' into ori_AlmeidaMEF
  file 'CayrouESC.peaks.mm10.bed'  into ori_CayrouESC

  script:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651111/suppl/GSM2651111%5FmES%5FWT%5Frepl%5FI%5Fpeaks%2EbigBed
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651112/suppl/GSM2651112%5FmES%5FWT%5Frepl%5FII%5Fpeaks%2EbigBed
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651107/suppl/GSM2651107%5FMEFs%5FWT%5Frepl%5FI%5Fpeaks%2EbigBed
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651108/suppl/GSM2651108%5FMEFs%5FWT%5Frepl%5FII%5Fpeaks%2EbigBed

  bigBedToBed GSM2651107_MEFs_WT_repl_I_peaks.bigBed   ori_MEF_rep1.gomez.tmp
  bigBedToBed GSM2651108_MEFs_WT_repl_II_peaks.bigBed  ori_MEF_rep2.gomez.tmp
  bigBedToBed GSM2651111_mES_WT_repl_I_peaks.bigBed    ori_ESC_rep1.gomez.tmp
  bigBedToBed GSM2651112_mES_WT_repl_II_peaks.bigBed   ori_ESC_rep2.gomez.tmp

  sort -k1,1 -k2n,2n -k3n,3n ori_MEF_rep1.gomez.tmp >ori_MEF_rep1.gomez.bed
  sort -k1,1 -k2n,2n -k3n,3n ori_MEF_rep2.gomez.tmp >ori_MEF_rep2.gomez.bed
  sort -k1,1 -k2n,2n -k3n,3n ori_ESC_rep1.gomez.tmp >ori_ESC_rep1.gomez.bed
  sort -k1,1 -k2n,2n -k3n,3n ori_ESC_rep2.gomez.tmp >ori_ESC_rep2.gomez.bed

  sort -k1,1 -k2n,2n -k3n,3n -m ori_MEF_rep1.gomez.bed ori_MEF_rep1.gomez.bed |mergeBed -i - >ori_MEF_gomez.bed
  sort -k1,1 -k2n,2n -k3n,3n -m ori_ESC_rep1.gomez.bed ori_ESC_rep1.gomez.bed |mergeBed -i - >ori_ESC_gomez.bed

  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651111/suppl/GSM2651111%5FmES%5FWT%5Frepl%5FI%2Ebw
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651112/suppl/GSM2651112%5FmES%5FWT%5Frepl%5FII%2Ebw
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651107/suppl/GSM2651107%5FMEFs%5FWT%5Frepl%5FI%2Ebw
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651108/suppl/GSM2651108%5FMEFs%5FWT%5Frepl%5FII%2Ebw

  bigWigToBedGraph GSM2651107_MEFs_WT_repl_I.bw  MEF_rep1.bedgraph
  bigWigToBedGraph GSM2651108_MEFs_WT_repl_II.bw MEF_rep2.bedgraph
  bigWigToBedGraph GSM2651111_mES_WT_repl_I.bw   ESC_rep1.bedgraph
  bigWigToBedGraph GSM2651112_mES_WT_repl_II.bw  ESC_rep2.bedgraph

  mapBed -a ori_MEF_gomez.bed -b MEF_rep2.bedgraph -c 4 -o sum >ori_MEF_gomez.bedgraph
  mapBed -a ori_ESC_gomez.bed -b ESC_rep2.bedgraph -c 4 -o sum >ori_ESC_gomez.bedgraph

  ## Peak calls
  wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651111/suppl/GSM2651111_mES_WT_repl_I_peaks.bigBed
  bigBedToBed GSM2651111_mES_WT_repl_I_peaks.bigBed GSM2651111_mES_WT_repl_I_peaks.bed
  sort -k1,1 -k2n,2n GSM2651111_mES_WT_repl_I_peaks.bed >AlmeidaESC.peaks.mm10.bed

  wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2651nnn/GSM2651107/suppl/GSM2651107_MEFs_WT_repl_I_peaks.bigBed
  bigBedToBed GSM2651107_MEFs_WT_repl_I_peaks.bigBed GSM2651107_MEFs_WT_repl_I_peaks.bed
  sort -k1,1 -k2n,2n GSM2651107_MEFs_WT_repl_I_peaks.bed >AlmeidaMEF.peaks.mm10.bed

  ## Cayrou_2015
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68347/suppl/GSE68347_Initiation_Sites.bedGraph.gz
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz

  gunzip GSE68347_Initiation_Sites.bedGraph.gz
  cut -f1-3  GSE68347_Initiation_Sites.bedGraph |sort -k1,1 -k2n,2n >CayrouESC.peaks.mm9.bed
  liftOver CayrouESC.peaks.mm9.bed  mm9ToMm10.over.chain.gz  CayrouESC.peaks.mm10.bed na

  """
  }

process getCpGIslands {

  publishDir params.outdirGomez, mode: 'copy', overwrite: true

  input:
  file (af) from af2.collect()
  file(winz) from winFiles_c.collect()

  output:
  file 'CGI.hg38.bedgraph'         into humanCGIBG
  file 'CGI.hg38.bed'              into humanCGIBED
  file 'CGI.mm10.bedgraph'         into mouseCGIBG
  file 'CGI.mm10.bed'              into mouseCGIBED
  file 'mm10_CpG.w1ks100.bedgraph' into (mm10CpGdensityBG_a, mm10CpGdensityBG_b)
  file 'mm10_CpG.w1ks100.bigwig'   into mm10CpGdensityBW
  file 'mm10_CpG_peaks.bed'        into (mm10CpGpeaksBED_a, mm10CpGpeaksBED_b)

  script:
  """
  wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz
  wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz

  wget http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-mm9.txt
  wget http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt

  cut -f1-3,6 model-based-cpg-islands-hg19.txt |grep -v start |sort -k1,1 -k2n,2n >CGI.hg19.bg
  cut -f1-3,6 model-based-cpg-islands-mm9.txt  |grep -v start |sort -k1,1 -k2n,2n >CGI.mm9.bg

  liftOver CGI.hg19.bg hg19ToHg38.over.chain.gz CGI.hg38.bedgraph na
  liftOver CGI.mm9.bg  mm9ToMm10.over.chain.gz  CGI.mm10.bedgraph na

  cut -f1-3 CGI.hg38.bedgraph >CGI.hg38.bed
  cut -f1-3 CGI.mm10.bedgraph >CGI.mm10.bed

  ## Get CpG density in 1 Kb wins 100bp slide
  perl ${params.codedir}/getDiNTsFromFA.pl ${mm10FA} CG >mm10_CpG.bed

  mapBed -a win1ks100.bed -b mm10_CpG.bed -g ${mm10IDX} -c 1 -o count |perl -lane 'print join("\\t",\$F[0],\$F[1]+450,\$F[2]-450,\$F[3])' >mm10_CpG.w1ks100.tmp
  sort -k1,1 -k2n,2n -k3n,3n mm10_CpG.w1ks100.tmp >mm10_CpG.w1ks100.bedgraph
  bedGraphToBigWig mm10_CpG.w1ks100.bedgraph ${mm10IDX} mm10_CpG.w1ks100.bigwig

  ## Get overlap stats of CpG peaks with Oris
  perl -lane '\$cpg = \$F[3]/1000; print \$_ if (\$cpg > 0.05)' mm10_CpG.w1ks100.bedgraph >CpG_thresholded.bed
  slopBed   -i CpG_thresholded.bed           -g ${mm10IDX} -b 100                                                   |sort -k1,1 -k2n,2n              >CpG_thresholded.slop100.bed
  mergeBed  -i CpG_thresholded.slop100.bed                                                                          |sort -k1,1 -k2n,2n              >CpG_thresholded.slopMerge.bed
  slopBed   -i CpG_thresholded.slopMerge.bed -g ${mm10IDX} -pct -r -0.5 -l -0.5 |slopBed -i - -g ${mm10IDX} -b 150  |sort -k1,1 -k2n,2n              >mm10_CpG_peaks.bed


  """
  }

process getAnnotationFiles {

  tag{sampleID}

  publishDir params.outdirAnnot, mode: 'copy', overwrite: true

  input:
  file (af) from af3.collect()

  output:
  file 'B6_maleHS.1bp.bedgraph'       into hotspotBG
  file 'B6_maleHS.500bp.bedgraph'     into hotspotBG500
  file 'B6_maleHS.oneMotif.500bp.bed' into hotspotOneMotif500a, hotspotOneMotif500b, hotspotOneMotif500c, hotspotOneMotif500d, hotspotOneMotif500e
  file 'B6_maleHS.1Kb.bed'            into hotspotBED
  file 'gencodeTSS.1bp.noMerge.bed'   into tssNoMergeBED1bp
  file 'gencodeTSS.1Kb.noMerge.bed'   into tssNoMergeBED1kb_a, tssNoMergeBED1kb_b
  file 'gencodeTSS.1KbDets.bed'       into tssBEDDets
  file 'gencodeTSS.1Kb.bed'           into tssBEDa  , tssBEDb  , tssBEDc  , tssBEDd  , tssBEDe  , tssBEDf  , tssBEDg  , tssBEDh, tssBEDi, tssBEDj, tssBEDk, tssBEDl, tssBEDm, tssBEDn
  file 'gencodeTES.1Kb.bed'           into gencodeTESBED
  file 'gencodeGene.bed'              into gencodeGeneBED
  file 'HS_and_TSS.1KbDets.bed'       into hstssBEDDets
  file 'HS_and_TSS.1Kb.bed'           into hstssBEDa, hstssBEDb, hstssBEDc, hstssBEDd, hstssBEDe, hstssBEDf, hstssBEDg, hstssBEDh, hstssBEDi, hstssBEDj, hstssBEDk
  file 'gencode.vM20.annotation.gtf'  into gencodeGTFa, gencodeGTFb, gencodeGTFc, gencodeGTFd
  file 'gencodeTSS.3Kb.noMerge.bed'   into tss3Kb
  file 'B6_maleHS.3Kb.bed'            into hs3Kb, hs3Kbb, hs3Kbc
  file 'B6_Pr9KO_maleHS.3Kb.bed'      into hsPrKO3Kb, hsPrKO3Kbb
  file 'refseqTSS.bed'                into refseqTSSBED,refseqTSSBEDa,refseqTSSBEDb
  file 'refseqTES.bed'                into refseqTESBED
  file 'refseqGene.bed'               into refseqGeneBED
  file 'refSeqEnhancers.bed'          into refseqEnhancersBED
  file 'repmaskerLINE.bed'            into rmLINEBED

  script:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664275/suppl/GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz
  gunzip -c GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_maleHS.bedgraph

  bedtools slop -l -0.5 -r -0.5 -pct -i B6_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6_maleHS.1bp.bedgraph
  cut -f1-3 B6_maleHS.1bp.bedgraph                                                                                       >B6_maleHS.1bp.bed

  bedtools slop -l 250  -r 250       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.500bp.bedgraph
  bedtools slop -l 500  -r 500       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.1Kb.bedgraph
  bedtools slop -l 500  -r 500       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6_maleHS.1Kb.bed
  bedtools slop -l 1500 -r 1500      -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6_maleHS.3Kb.bed
  cat B6_maleHS.1Kb.bedgraph |perl -lane \'print join("\\t",@F[0..3],"HS","+")\'                 >B6_maleHS.1Kb.forMerge.bed

  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664291/suppl/GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz
  gunzip -c GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_Pr9KO_maleHS.bedgraph

  bedtools slop -l -0.5 -r -0.5 -pct -i B6_Pr9KO_maleHS.bedgraph     -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6_Pr9KO_maleHS.1bp.bedgraph
  cut -f1-3 B6_Pr9KO_maleHS.1bp.bedgraph                                                                                           >B6_Pr9KO_maleHS.1bp.bed

  bedtools slop -l 250  -r 250       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX}            >B6_Pr9KO_maleHS.500bp.bedgraph
  bedtools slop -l 500  -r 500       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX}            >B6_Pr9KO_maleHS.1Kb.bedgraph
  bedtools slop -l 500  -r 500       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX} |cut -f1-3 >B6_Pr9KO_maleHS.1Kb.bed
  bedtools slop -l 1500 -r 1500      -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX} |cut -f1-3 >B6_Pr9KO_maleHS.3Kb.bed

  perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' B6_maleHS.500bp.bedgraph >B6HS500forFIMO.bed
  bedtools getfasta -fi ${mm10FA} -bed B6HS500forFIMO.bed -name -fo B6_maleHS.500bp.fa

  fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo ${params.datadir}/PRDM9motif/PRBS_B6.MEMEv4.pwm B6_maleHS.500bp.fa

  perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo ./fimo/fimo.tsv --w 250 --out B6_maleHS.oneMotif.500bp.bed

  ##GENCODE
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz
  ##TSS
  perl ${params.codedir}/gencodeGTFtoTSS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >gencodeTSS.1bp.noMerge.bed

  bedtools slop -l 500 -r 500   -i gencodeTSS.1bp.noMerge.bed -g ${mm10IDX} >gencodeTSS.1Kb.noMerge.bed
  bedtools slop -l 1500 -r 1500 -i gencodeTSS.1bp.noMerge.bed -g ${mm10IDX} >gencodeTSS.3Kb.noMerge.bed
  mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX}            >gencodeTSS.1KbDets.bed
  mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |cut -f1-3 >gencodeTSS.1Kb.bed

  cat gencodeTSS.1KbDets.bed |perl -lane \'print join("\\t",@F[0..2],0,@F[4..5])\'                           >gencodeTSS.1Kb.forMerge.bed

  ##TES
  perl ${params.codedir}/gencodeGTFtoTES.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |grep -P "\\s+" >gencodeTES.1bp.noMerge.bed

  bedtools slop -l 500 -r 500   -i gencodeTES.1bp.noMerge.bed -g ${mm10IDX} >gencodeTES.1Kb.noMerge.bed
  mergeBed -i gencodeTES.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |\
                                            ${params.codedir}/sortBEDByFAI.pl - \
                                            ${mm10IDX} >gencodeTES.1Kb.bed

  ##GENES
  perl ${params.codedir}/gencodeGTFtoCDS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >gencodeGene.noMerge.bed
  mergeBed -s -i gencodeGene.noMerge.bed -c 4,5,6 \
  -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - \
            ${mm10IDX} |grep -P "\\s+" >gencodeGene.bed

  ##HS and TSS
  cat B6_maleHS.1Kb.forMerge.bed gencodeTSS.1Kb.forMerge.bed |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX}            >HS_and_TSS.1KbDets.bed
  cat B6_maleHS.1Kb.forMerge.bed gencodeTSS.1Kb.forMerge.bed |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |cut -f1-3 >HS_and_TSS.1Kb.bed

  gunzip gencode.vM20.annotation.gtf.gz

  ##REFSEQ
  wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqCurated.txt.gz
  perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz TSS  |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >refseqTSS.bed
  perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz TES  |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >refseqTES.bed
  perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz Gene |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >refseqGene.bed

  ##ENHANCERS
  bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/mm10/ncbiRefSeq/refSeqFuncElems.bb stdout |cut -f1-3 |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >refSeqEnhancers.bed

  ##LINES
  wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
  zcat rmsk.txt.gz | perl -lane 'print join("\\t",@F[5..7],\$F[10],\$F[11],\$F[9]) if (\$F[1] > 20000 && \$_ =~ /LINE/)' |grep -P "\\s+" |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >repmaskerLINE.bed

  """
  }

process convertReplicationDomainsOrgBGs_Mouse {

  publishDir params.outdirBG, mode: 'copy', overwrite: true

  input:
  file(af4)
  file(bg) from rdOrgBG_mm10.collect()

  output:
  file '*forModel.bedgraph' into rdOrgBG_forModel_mm10

  script:
  """
  bedtools makewindows -g ${mm10IDX} -w 500000 -s 50000  | perl -lane 'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)' |sort -k1,1 -k2n,2n -k3n,3n >win500k50k.bed

  for bg in *bedgraph; do

    nm=\${bg/.bedgraph/}
    modBG=\${bg/bedgraph/forModel.bedgraph}

    grep -P '^chr' \$bg |sort -k1,1 -k2n,2n -k3rn,3rn >\$nm.s1.bg

    mapBed -a win500k50k.bed -b \$nm.s1.bg -sorted -c 4,4 -o mean,count |\
        perl -lane '\$F[3] = ((\$F[3] eq ".")?"0":\$F[3]); print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[2]-250000+24999,\$F[3]) if (\$F[4] > 30)' >\$nm.s2.bg

    perl ${params.codedir}/convertToModellingBG.pl \$nm.s2.bg ${mm10IDX} >\$nm.s3.bg

    intersectBed -a \$nm.s3.bg -b ${params.datadir}/blacklist/mm10.blacklist.bed  -v >\$modBG
  done
  """
  }

process convertReplicationDomainsOrgBGs_Human {

  publishDir params.outdirBG, mode: 'copy', overwrite: true

  input:
  file (af) from af5.collect()
  file(bg) from rdOrgBG_hg38.collect()

  output:
  file '*forModel.bedgraph' into rdOrgBG_forModel_hg38

  script:
  """
  wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

  wget -O RPE1_MID1_RT.hg19.bg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2904948&format=file&file=GSM2904948%5FP285%5F06%5F1%5Fw200ks40k%5Fmap%5Fcount%5Fmedian%5Flog2%2EbedGraph%2Egz"
  wget -O RPE1_MID2_RT.hg19.bg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2904949&format=file&file=GSM2904949%5FP285%5F07%5F1%5Fw200ks40k%5Fmap%5Fcount%5Fmedian%5Flog2%2EbedGraph%2Egz"
  wget -O RPE1_MID3_RT.hg19.bg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2904950&format=file&file=GSM2904950%5FP285%5F08%5F1%5Fw200ks40k%5Fmap%5Fcount%5Fmedian%5Flog2%2EbedGraph%2Egz"

  gunzip RP*gz

  liftOver RPE1_MID1_RT.hg19.bg  hg19ToHg38.over.chain.gz RPE1_MID1_RT_hg38.bedgraph na
  liftOver RPE1_MID2_RT.hg19.bg  hg19ToHg38.over.chain.gz RPE1_MID2_RT_hg38.bedgraph na
  liftOver RPE1_MID3_RT.hg19.bg  hg19ToHg38.over.chain.gz RPE1_MID3_RT_hg38.bedgraph na

  bedtools makewindows -g ${hg38IDX} -w 500000 -s 50000  | perl -lane 'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)' |sort -k1,1 -k2n,2n -k3n,3n >win500k50k.bed
  for bg in *bedgraph; do

    nm=\${bg/.bedgraph/}
    modBG=\${bg/bedgraph/forModel.bedgraph}

    grep -P '^chr' \$bg |sort -k1,1 -k2n,2n -k3rn,3rn >\$nm.s1.bg

    pat="RPE1"

    if [[ \$nm =~ \$pat ]]; then
      mapBed -a win500k50k.bed -b \$nm.s1.bg -sorted -c 4,4 -o mean,count |\
      perl -lane '\$F[3] = ((\$F[3] eq ".")?"0":\$F[3]); print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[2]-250000+24999,\$F[3]) if (\$F[4] > 9)' >\$nm.s2.bg
    else
      mapBed -a win500k50k.bed -b \$nm.s1.bg -sorted -c 4,4 -o mean,count |\
      perl -lane '\$F[3] = ((\$F[3] eq ".")?"0":\$F[3]); print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[2]-250000+24999,\$F[3]) if (\$F[4] > 30)' >\$nm.s2.bg
    fi

    perl ${params.codedir}/convertToModellingBG.pl \$nm.s2.bg ${hg38IDX} >\$nm.s3.bg

    intersectBed -a \$nm.s3.bg -b ${params.datadir}/blacklist/hg38.blacklist.bed  -v >\$modBG
  done
  """
  }

// process getOrigins {
//
//   publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*svg'
//   publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*png'
//   publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*pdf'
//   publishDir params.outdirOrigins, mode: 'copy', overwrite: true, pattern: '*bedgraph'
//   publishDir params.outdirOrigins, mode: 'copy', overwrite: true, pattern: '*bed'
//   publishDir params.outdirOrigins, mode: 'copy', overwrite: true, pattern: '*tab'
//   publishDir params.outdirOrigins, mode: 'copy', overwrite: true, pattern: '*Rdata'
//
//   input:
//   file (af) from af6.collect()
//   //set val(ssdsFolder), val(genome), val(outName) from oriSSDSPaths
//
//   output:
//   file '*mm10_OriSSDS*.bed'           into mouseOriginsInitBED
//   file '*mm10_OriSSDS*.tab'           into mouseOriginsInitTAB
//   file '*png'                         into oriPNGFigs
//   file '*pdf'                         into oriPDFFigs
//   file "hiconf_origins.mm10.bedgraph" into (mouseOriginsBG_a, mouseOriginsBG_b, mouseOriginsBG_c, mouseOriginsBG_d, mouseOriginsBG_e)
//   file "rand*hiconf_*.mm10.bedgraph"  into mouseRandomOriginsBG
//   file "union_origins.mm10.bedgraph"  into mouseOriginsBGUnion
//   file "WT4_origins.mm10.bedgraph"    into mouseOriginsWT4BG
//   file '*.Rdata'                      into (allOriginRData,allOriginRData_a)
//   file '*rigins*bed*'                 into (allOriginsBED_a, allOriginsBED_b)
//
//
//   script:
//   """
//   ## Get Project .Rprofile file
//   cp accessoryFiles/scripts/R/Rprofile.workflow    ./.Rprofile
//   echo RTSCRIPTS="accessoryFiles/scripts/R/"     >>./.Renviron
//
//   profiles=`echo "${workflow.profile}" |perl -pi -e 's/modules/modules,callOrigins/'`;
//
//   env NXF_GENOMES="accessoryFiles/genomeFiles/" nextflow run ${params.codedir}/callOriginsSNS.nf \
//                                                        -c accessoryFiles/config/nextflow.config.nf \
//                                                        -profile \$profiles \
//                                                        --projectdir `pwd` \
//                                                        --outdir . \
//                                                        --ssdsPath accessoryFiles/data/oriSSDS/mouse/ \
//                                                        --genome mm10 \
//                                                        --name ${outName} \
//                                                        --extSize 2000 \
//                                                        --reps 1 \
//                                                        --blacklist accessoryFiles/data/blacklist/mm10.blacklist.bed \
//                                                        --accessoryFiles accessoryFiles
//
//   if [[ ${genome} == "mm10" ]]; then
//     R --no-save <accessoryFiles/scripts/R/distillOriginsMM10.R ||true
//     bedtools shuffle -chrom -i hiconf_origins.mm10.bedgraph -g ${mm10IDX} |sort -k1,1 -k2n,2n -k3n,3n >randomized_hiconf_origins.mm10.bedgraph
//   fi
//
//   """
//   }

/// Calling Origins ////////////////////////////////////////////////////////////
process shufBEDs {

  tag { ssdsBed }

  input:
  set val(id), val(genome), file(ssdsBed) from ssdsBEDs_a
  file (af)                               from af31.collect()

  output:
  file "*.T.sq30.bed" into sqtBed, sqtBed2, sqtBed3, sqtBed4, sqtBed5

  script:
  """
  echo ${ssdsBed}
  perl accessoryFiles/scripts/q30SSDS.pl <${ssdsBed} >${id}.T.q30.bed
  intersectBed -a ${id}.T.q30.bed -b  accessoryFiles/data/blacklist/mm10.blacklist.bed -v |grep -v chrM >${id}.T.q30noBL.bed
  shuf ${id}.T.q30noBL.bed |perl accessoryFiles/scripts/reverseStrandsForOriCalling.pl >${id}.T.sq30.bed
  """
  }

process callReplicationOrigins {

   tag { shufPC }

   input:
   file (sqtBed)
   file (af) from af32.collect()
   each shufPC from satCurvePCs

   output:
   file "*.peaks.bed" into satCurvePeaks,scPeaks
   file "*.peaks.xls" into scXLS

   script:
   def sqtNm = sqtBed.name.replaceFirst(".T.sq30.bed","")
   """
   nT=`cat ${sqtBed} |wc -l`
   nPC=`perl -e 'print int('\$nT'*${shufPC})'`

   perl accessoryFiles/scripts/pickNlines.pl ${sqtBed} \$nPC >\$nPC.tmp
   sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq >\$nPC.T.bed

   touch ${sqtNm}'.N'\$nPC'_'${shufPC}'pc.999_peaks.xls'

   for i in {0..${params.reps}}; do
     thisName=${sqtNm}'.N'\$nPC'_${shufPC}pc.'\$i
     perl accessoryFiles/scripts/pickNlines.pl ${sqtBed} \$nPC >\$nPC.tmp

     sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 \$nPC.tmp |uniq >\$nPC.T.bed

     macs2 callpeak  -g mm \\
       -t \$nPC.T.bed \\
       -q 0.001 \\
       --extsize 2000 \\
       --nomodel  \\
       --nolambda \\
       --name ${sqtNm}

     cut -f1-3 ${sqtNm}_peaks.narrowPeak |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >\$thisName.peaks.bed
     mv ${sqtNm}_peaks.xls \$thisName.peaks.xls
   done
   rm *999_peaks.xls
   """
  }

process makeOriginCallingSaturationCurve {

   publishDir params.outdirOrigins,  mode: 'copy', overwrite: true, pattern: '*png'

   input:
   file (af) from af33.collect()
   file (sqtBed2)
   file (pk) from satCurvePeaks.collect()

   output:
   file "*satCurve.tab"                           into satCurveTable
   file "*.png"                                   into satCurvePNG
   file "*preditedHotspotsFromMoreSequencing.tab" into preditedHS

   script:
   def sqtNm = sqtBed2.name.replaceFirst(".T.sq30.bed","")
   """
   ## Get Project .Rprofile file
   cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

   wc -l ${sqtNm}*peaks.bed |grep -v total >${sqtNm}_allVals.tab
   perl accessoryFiles/scripts/drawMACSSaturationPlot.pl ${sqtNm}_allVals.tab ${sqtNm}.satCurve.tab ${sqtNm}
   rm ${sqtNm}_allVals.tab
   """

  }

process processAndKeepOrigins {

   publishDir params.outdirOrigins,  mode: 'copy', overwrite: true, pattern: "*peaks.RC.bedgraph"

   input:
   file (af) from af34.collect()
   file (pk)      from scPeaks.collect()
   file (pkXLS)   from scXLS.collect()
   file (sqtBed3)

   output:
   file "*peaks.RC.bed"      into finalPeaks, fPeaks2
   file "*peaks.RC.bedgraph" into finalBG
   file "*finalpeaks.xls"    into finalXLS

   script:
   def sqtNm = sqtBed3.name.replaceFirst(".T.sq30.bed","")
   """
   sort -k1,1 -k2n,2n -k3n,3n ${sqtNm}*_1.00pc.*.peaks.bed >allpks.bed
   mergeBed -i allpks.bed -d 500 >${sqtNm}.peaks.bedM
   perl accessoryFiles/scripts/normalizeStrengthByAdjacentRegions.pl --bed ${sqtNm}.peaks.bedM --in ${sqtBed3} --out ${sqtNm}.peaks.RC.bedgraph --rc --tmp tmp
   cut -f1-3 ${sqtNm}.peaks.RC.bedgraph >${sqtNm}.peaks.RC.bed
   cat ${sqtNm}*1.00pc.0.peaks.xls >${sqtNm}.finalpeaks.xls
   """
  }

process mergeOrigins {

   //publishDir params.outdir,  mode: 'copy', overwrite: true

   input:
   file(pk)   from finalPeaks.collect()
   file (af) from af35.collect()

   output:
   file "*mergedpeaks.tab" into (finalMergedPeaks_a,finalMergedPeaks_b,mergeA)
   file "aP.nuc"           into mergednuc

   script:
   """
   sort -k1,1 -k2n,2n -k3n,3n *peaks.RC.bed >allPeaks.set
   mergeBed -i allPeaks.set -d -500 >all.mergedpeaks.tab
   cut -f1-3 all.mergedpeaks.tab >aP.bed

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GC -C -bed aP.bed |cut -f4-5,13 |grep -v pct >nGC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GA -C -bed aP.bed |cut -f13  |grep -v user >nGA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GG -C -bed aP.bed |cut -f13  |grep -v user >nGG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GT -C -bed aP.bed |cut -f13  |grep -v user >nGT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CC -C -bed aP.bed |cut -f13  |grep -v user >nCC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CA -C -bed aP.bed |cut -f13  |grep -v user >nCA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CG -C -bed aP.bed |cut -f13  |grep -v user >nCG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CT -C -bed aP.bed |cut -f13  |grep -v user >nCT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AC -C -bed aP.bed |cut -f13  |grep -v user >nAC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AA -C -bed aP.bed |cut -f13  |grep -v user >nAA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AG -C -bed aP.bed |cut -f13  |grep -v user >nAG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AT -C -bed aP.bed |cut -f13  |grep -v user >nAT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TC -C -bed aP.bed |cut -f13  |grep -v user >nTC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TA -C -bed aP.bed |cut -f13  |grep -v user >nTA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TG -C -bed aP.bed |cut -f13  |grep -v user >nTG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TT -C -bed aP.bed |cut -f13  |grep -v user >nTT

   echo -e 'AT GC GpC GpG GpA GpT CpC CpG CpA CpT ApC ApA ApG ApT TpC TpA TpG TpT' |perl -pi -e \'s/\\s+(\\S)/\\t\$1/g\' >aP.ncALL
   paste nGC nGG nGA nGT nCC nCG nCA nCT nAC nAA nAG nAT nTC nTA nTG nTT  >>aP.ncALL

   bedtools slop -l -0.5 -r -0.5 -i aP.bed  -pct -g accessoryFiles/genomeFiles//mm10/genome.fa.fai >aPcenter.bed
   bedtools slop -l 250  -r 250  -i aPcenter.bed   -g accessoryFiles/genomeFiles//mm10/genome.fa.fai >w500.bed

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GC -C -bed w500.bed |cut -f4-5,13 |grep -v pct >wGC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GA -C -bed w500.bed |cut -f13  |grep -v user >wGA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GG -C -bed w500.bed |cut -f13  |grep -v user >wGG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern GT -C -bed w500.bed |cut -f13  |grep -v user >wGT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CC -C -bed w500.bed |cut -f13  |grep -v user >wCC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CA -C -bed w500.bed |cut -f13  |grep -v user >wCA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CG -C -bed w500.bed |cut -f13  |grep -v user >wCG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern CT -C -bed w500.bed |cut -f13  |grep -v user >wCT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AC -C -bed w500.bed |cut -f13  |grep -v user >wAC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AA -C -bed w500.bed |cut -f13  |grep -v user >wAA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AG -C -bed w500.bed |cut -f13  |grep -v user >wAG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern AT -C -bed w500.bed |cut -f13  |grep -v user >wAT

   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TC -C -bed w500.bed |cut -f13  |grep -v user >wTC
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TA -C -bed w500.bed |cut -f13  |grep -v user >wTA
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TG -C -bed w500.bed |cut -f13  |grep -v user >wTG
   nucBed -fi accessoryFiles/genomeFiles/mm10/genome.fa -pattern TT -C -bed w500.bed |cut -f13  |grep -v user >wTT

   echo -e 'AT500 GC500 GpC500 GpG500 GpA500 GpT500 CpC500 CpG500 CpA500 CpT500 ApC500 ApA500 ApG500 ApT500 TpC500 TpA500 TpG500 TpT500' |perl -pi -e \'s/\\s+(\\S)/\\t\$1/g\' >aP.nc500
   paste wGC wGG wGA wGT wCC wCG wCA wCT wAC wAA wAG wAT wTC wTA wTG wTT >>aP.nc500

   paste aP.ncALL aP.nc500 >aP.nuc
   """
  }

process getOriginOverlaps {

   //publishDir params.outdir,  mode: 'copy', overwrite: true

   input:
   file (fPeaks2)
   file (mergeA)

   output:
   file "*.OL" into mergedOL

   script:
   def nm = fPeaks2.name.replaceFirst(".peaks.RC.bed","")
   """

   echo "ori_"${nm} >${fPeaks2}.OL
   intersectBed -a ${mergeA} -b ${fPeaks2} -c |cut -f4 >>${fPeaks2}.OL
   """
  }

process recalcOriginStrength {

   //publishDir params.outdir,  mode: 'copy', overwrite: true

   input:
   file (af) from af36.collect()
   file (sqtBed5)
   file (mBed)    from finalMergedPeaks_a

   output:
   file "*finalmergedpeaks.col" into mergedStr
   file "*finalmergedpeaks.FFRR" into mergedFFRR

   script:
   def sqtNm = sqtBed5.name.replaceFirst(".T.sq30.bed","")
   """
   perl accessoryFiles/scripts/normalizeStrengthByAdjacentRegions.pl --bed ${mBed} --in ${sqtBed5} --out ${sqtNm}.finalmergedpeaks.initout --tmp tmp
   cut -f4 ${sqtNm}.finalmergedpeaks.tab |perl -pi -e \'s/strength/str_\'${sqtNm}\'/\' >${sqtNm}.finalmergedpeaks.col
   cut -f9-12 ${sqtNm}.finalmergedpeaks.tab  >${sqtNm}.finalmergedpeaks.FFRR
   perl -pi -e \'s/leftFwd/LW_\'${sqtNm}\'/\'  ${sqtNm}.finalmergedpeaks.FFRR
   perl -pi -e \'s/rightFwd/RW_\'${sqtNm}\'/\' ${sqtNm}.finalmergedpeaks.FFRR
   perl -pi -e \'s/leftRev/LC_\'${sqtNm}\'/\'  ${sqtNm}.finalmergedpeaks.FFRR
   perl -pi -e \'s/rightRev/RC_\'${sqtNm}\'/\' ${sqtNm}.finalmergedpeaks.FFRR
   """
  }

process finalMergeForOrigins {

   publishDir params.outdirOrigins,  mode: 'copy', overwrite: true

   input:
   file (af) from af37.collect()
   file(mBed)   from finalMergedPeaks_b
   file(sCol)   from mergedStr.collect()
   file(frCol)  from mergedFFRR.collect()
   file(mOL)    from mergedOL.collect()
   file(mNuc)   from mergednuc.collect()

   output:
   file '*mm10_OriSSDS*.bed'           into mouseOriginsInitBED
   file '*mm10_OriSSDS*.tab'           into mouseOriginsInitTAB
   file "hiconf_origins.mm10.bedgraph" into (mouseOriginsBG_a, mouseOriginsBG_b, mouseOriginsBG_c, mouseOriginsBG_d, mouseOriginsBG_e, mouseOriginsBG_f, mouseOriginsBG_g, mouseOriginsBG_h)
   file "hiconf_origins.mm10.Rdata"    into mouseOriginsHCRdata
   file "rand*hiconf_*.mm10.bedgraph"  into mouseRandomOriginsBG
   file "union_origins.mm10.bedgraph"  into mouseOriginsBGUnion
   file "WT4_origins.mm10.bedgraph"    into mouseOriginsWT4BG
   file '*.Rdata'                      into (allOriginRData,allOriginRData_a)
   file '*rigins*bed*'                 into (allOriginsBED_a, allOriginsBED_b)

   script:
   def sqtName = "mm10_OriSSDS"
   """
   echo -e \'cs\\tfrom\\tto\' >m.head
   cat ${mBed} >>m.head
   paste m.head *OL *col *FFRR *nuc >${sqtName}.Origins.tmp
   perl -lane \'if (\$_ =~ /from/){for my \$i(0..\$#F){if (\$i > 0){\$F[\$i] = \$F[\$i].(\$Flst{\$F[\$i]}++?"_".\$Flst{\$F[\$i]}:"")}}}; print join("\\t",@F)\' ${sqtName}.Origins.tmp >${sqtName}.Origins.tab
   cut -f1-3 ${sqtName}.Origins.tab |grep -v from >${sqtName}.Origins.bed

   R --no-save <accessoryFiles/scripts/R/distillOriginsMM10.R ||true
   bedtools shuffle -chrom -i hiconf_origins.mm10.bedgraph -g ${mm10IDX} |sort -k1,1 -k2n,2n -k3n,3n >randomized_hiconf_origins.mm10.bedgraph
   """
  }
/// END: Calling Origins ///////////////////////////////////////////////////////
//
// process getBendabilityAtOrigins {
//
//   publishDir params.outdirFigs  , mode: 'copy', overwrite: true, pattern: '*png'
//   publishDir params.outdirFigs  , mode: 'copy', overwrite: true, pattern: '*pdf'
//
//   input:
//   file (af) from af7.collect()
//   file(oriMM) from mouseOriginsBG_d
//
//   output:
//   file "*png"              into oriBendabilityPNG
//   file "*pdf"              into oriBendabilityPDF
//   file "*.bendability.txt" into oriBendabilityTxt
//
//   script:
//   """
//   cp ${params.codedir}/CalcPhysicoChemicalPropsFromFA.pl .
//   cp ${params.codeRdir}/plotBendability.R .
//
//   bedtools slop -l -0.5 -r -0.5 -pct -i ${oriMM} -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |cut -f1-3 >oriMMnt.bed
//
//   bedtools slop -l 10000 -r 10000 -i oriMMnt.bed -g ${mm10IDX} >oriMM.10k.bed
//
//   bedtools getfasta -fi ${mm10FA} -bed oriMM.10k.bed -fo oriMM.10k.fa
//
//   perl CalcPhysicoChemicalPropsFromFA.pl oriMM.10k.fa 25 >mm10Origins.bendability.txt
//
//   Rscript plotBendability.R -i mm10Origins.bendability.txt -n mm10_hiconf_origins
//   """
//   }
//
process getCoverageAtOrigins {

  tag  {bed}

  publishDir params.outdirBW    , mode: 'copy', overwrite: true, pattern: '*bigwig'
  publishDir params.outdirDTdata, mode: 'copy', overwrite: true, pattern: '*gz'

  input:
  file (af) from af8.collect()
  set val(id), val(genome), file(bed) from ssdsBEDs
  file(ori)                           from allOriginsBED_a.collect()
  file(winz)                          from winFiles_e

  output:
  file "*RPKM.bigwig"             into (oriSSDSBW1, oriSSDSBW2)
  file "*FR.frags.bigwig"         into (oriSSDSFRbw)
  file "*.origins*fragsMatrix.gz" into (oriSSDS_dtMatrix,oriSSDS_dtMatrix_a,oriSSDS_dtMatrix_b)
  file "*.g4.*.fragsMatrix.gz"    into g4SSDS_dtMatrix
  file "*png"                     into g4heatmapFigs

  script:
  """
  iName="${genome}_${id}"
  oName=\${iName/.bed/}

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bed} |uniq >noDuplicates.bed
  grep -P  \'\\+\' noDuplicates.bed >\$oName.POS.frags.bed
  grep -vP \'\\+\' noDuplicates.bed >\$oName.NEG.frags.bed

  bedtools bedtobam -i \$oName.POS.frags.bed -g ${params.genomedir}/${genome}_genome.fa.fai >\$oName.POS.frags.initbam
  samtools sort --reference ${params.genomedir}/${genome}_genome.fa.fai \$oName.POS.frags.initbam >\$oName.POS.frags.bam
  samtools index \$oName.POS.frags.bam

  bedtools bedtobam -i \$oName.NEG.frags.bed -g ${params.genomedir}/${genome}_genome.fa.fai >\$oName.NEG.frags.initbam
  samtools sort --reference ${params.genomedir}/${genome}_genome.fa.fai \$oName.NEG.frags.initbam >\$oName.NEG.frags.bam
  samtools index \$oName.NEG.frags.bam

  bedtools bedtobam -i noDuplicates.bed -g ${params.genomedir}/${genome}_genome.fa.fai >\$oName.frags.initbam
  samtools sort --reference ${params.genomedir}/${genome}_genome.fa.fai \$oName.frags.initbam >\$oName.ALL.frags.bam
  samtools index \$oName.ALL.frags.bam

  bamCoverage --bam \$oName.POS.frags.bam --outFileName \$oName.POS.frags.RPKM.bigwig --outFileFormat bigwig --binSize 150 -p max --ignoreForNormalization chrX chrM chrY --normalizeUsing RPKM
  bamCoverage --bam \$oName.NEG.frags.bam --outFileName \$oName.NEG.frags.RPKM.bigwig --outFileFormat bigwig --binSize 150 -p max --ignoreForNormalization chrX chrM chrY --normalizeUsing RPKM
  bamCoverage --bam \$oName.ALL.frags.bam --outFileName \$oName.ALL.frags.RPKM.bigwig --outFileFormat bigwig --binSize 150 -p max --ignoreForNormalization chrX chrM chrY --normalizeUsing RPKM

  mapBed -a win1ks100.bed -b \$oName.POS.frags.bam -c 1 -o count >\$oName.POS.frags.bg
  mapBed -a win1ks100.bed -b \$oName.NEG.frags.bam -c 1 -o count >\$oName.NEG.frags.bg

  paste \$oName.POS.frags.bg \$oName.NEG.frags.bg |perl -lane 'use Math::Round; \$fr = (\$F[3]+1)/(\$F[7]+1); \$logfr = log(\$fr)/log(2); \$mid=round((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$mid-49,\$mid+50,\$logfr)' |sort -k1,1 -k2n,2n >\$oName.FR.frags.bg
  bedGraphToBigWig \$oName.FR.frags.bg ${params.genomedir}/${genome}_genome.fa.fai \$oName.FR.frags.bigwig

  computeMatrix reference-point --referencePoint center -a 10000 -b 10000 \
                        -R ${genome}_OriSSDS.Origins.bed \
                        -S \$oName.POS.frags.RPKM.bigwig \
                        \$oName.NEG.frags.RPKM.bigwig \
                        \$oName.ALL.frags.RPKM.bigwig \
                        -p max -o \$oName.origins.fragsMatrix.gz \
                        --binSize 100

  intersectBed -a ${params.datadir}/g4/${genome}.qparser2.g4.bed -b ${genome}_OriSSDS.Origins.bed -v |sort -k1,1 -k2n,2n >g4noOri.bed
  perl -lane \'@X=split(":",\$F[4]); print join("\\t",\$F[0],\$F[1]-1,\$F[1]) if (\$X[0] >= 6 && \$F[5] eq "+")\' g4noOri.bed >g4.Watson.bed
  perl -lane \'@X=split(":",\$F[4]); print join("\\t",\$F[0],\$F[2]-1,\$F[2]) if (\$X[0] >= 6 && \$F[5] ne "+")\' g4noOri.bed >g4.Crick.bed

  computeMatrix reference-point --referencePoint center -a 10000 -b 10000 \
                        -R g4.Watson.bed \
                        -S \$oName.POS.frags.RPKM.bigwig \
                        \$oName.NEG.frags.RPKM.bigwig \
                        \$oName.ALL.frags.RPKM.bigwig \
                        -p max -o \$oName.g4.Watson.fragsMatrix.gz \
                        --binSize 100

  computeMatrix reference-point --referencePoint center -a 10000 -b 10000 \
                        -R g4.Crick.bed \
                        -S \$oName.POS.frags.RPKM.bigwig \
                        \$oName.NEG.frags.RPKM.bigwig \
                        \$oName.ALL.frags.RPKM.bigwig \
                        -p max -o \$oName.g4.Crick.fragsMatrix.gz \
                        --binSize 100

  bamCoverage --bam \$oName.POS.frags.bam --outFileName \$oName.POS.frags.25.RPKM.bigwig --outFileFormat bigwig --binSize 25 -p max --ignoreForNormalization chrX chrM chrY
  bamCoverage --bam \$oName.NEG.frags.bam --outFileName \$oName.NEG.frags.25.RPKM.bigwig --outFileFormat bigwig --binSize 25 -p max --ignoreForNormalization chrX chrM chrY

  computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
  -R g4.Watson.bed g4.Crick.bed \
  -S \$oName.POS.frags.25.RPKM.bigwig \
  -p max -o \$oName.g4x2Watson.25.fragsMatrix.gz \
                                --binSize 100

  computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
  -R g4.Watson.bed g4.Crick.bed \
  -S \$oName.NEG.frags.25.RPKM.bigwig \
  -p max -o \$oName.g4x2Crick.25.fragsMatrix.gz \
                                --binSize 100

  plotHeatmap --matrixFile \$oName.g4x2Watson.25.fragsMatrix.gz \
  -o \$oName.g4x2Watson.png --xAxisLabel "Distance to g4 (Kb)" \
  --heatmapWidth 5 --colorMap Reds --samplesLabel "Watson" \
  --dpi 300 --heatmapHeight 10 \
              --regionsLabel "Top-strand g4s" "Bottom-strand g4s" \
              --refPointLabel g4 --averageTypeSummaryPlot  median \
              --whatToShow "heatmap only"

  plotHeatmap --matrixFile \$oName.g4x2Crick.25.fragsMatrix.gz \
  -o \$oName.g4x2Crick.png --xAxisLabel "Distance to g4 (Kb)" \
  --heatmapWidth 5 --colorMap Reds --samplesLabel "Crick" \
  --dpi 300 --heatmapHeight 10 \
              --regionsLabel "Top-strand g4s" "Bottom-strand g4s" \
              --refPointLabel g4 --averageTypeSummaryPlot  median \
              --whatToShow "heatmap only"

  if [ ${genome} == 'mm10' ]; then
    cut -f1-3 hiconf_origins.mm10.bedgraph >hiconf_origins.mm10.bed
    computeMatrix reference-point --referencePoint center -a 10000 -b 10000 \
                        -R hiconf_origins.mm10.bed \
                        -S \$oName.POS.frags.RPKM.bigwig \
                        \$oName.NEG.frags.RPKM.bigwig \
                        \$oName.ALL.frags.RPKM.bigwig \
                        -p max -o \$oName.hiconfOri.mm10.fragsMatrix.gz \
                        --binSize 100

    cut -f1-3 WT4_origins.mm10.bedgraph >WT4_origins.mm10.bed
    computeMatrix reference-point --referencePoint center -a 10000 -b 10000 \
                        -R WT4_origins.mm10.bed \
                        -S \$oName.POS.frags.RPKM.bigwig \
                        \$oName.NEG.frags.RPKM.bigwig \
                        \$oName.ALL.frags.RPKM.bigwig \
                        -p max -o \$oName.WT4Ori.mm10.fragsMatrix.gz \
                        --binSize 100
  fi

  """
  }

process plotCoverageVG4 {

  publishDir params.outdirFigs, mode: 'copy', overwrite: true

  input:
  file (af) from af9.collect()
  file(oriGZ) from oriSSDS_dtMatrix.collect()
  file(g4GZ)  from g4SSDS_dtMatrix.collect()

  output:
  file '*png'                        into oriVG4PNGFigs
  file '*pdf'                        into oriVG4PDFFigs

  script:
  """
  for gz in *.fragsMatrix.gz; do
    gout=\${gz/.gz/}
    gunzip -c \$gz >\$gout
  done

  genome="mm10"

  intersectBed -a ${params.datadir}/g4/\$genome".qparser2.g4.bed" -b \$genome"_OriSSDS.Origins.bed" -v |sort -k1,1 -k2n,2n >g4noOri.bed
  perl -lane \'@X=split(":",\$F[4]); print join("\\t",\$F[0],\$F[1]-1,\$F[1]) if (\$X[0] >= 6 && \$F[5] eq "+")\' g4noOri.bed >g4.Watson.bed
  perl -lane \'@X=split(":",\$F[4]); print join("\\t",\$F[0],\$F[2]-1,\$F[2]) if (\$X[0] >= 6 && \$F[5] ne "+")\' g4noOri.bed >g4.Crick.bed

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  R --no-save <accessoryFiles/scripts/R/plot_Origins_v_G4_Coverage.R ||true
  """
  //R --no-save --no-site-file --no-init-file --no-restore --silent --slave <${params.codeRdir}/plot_Origins_v_G4_Coverage_hg38.R ||true

  }

// MAKE SLICE DATA
Channel.from( ["chr1" , 38764000,  39764000, "chr1_38Mb_figure1"],
              ["chr2" , 74694564,  74768877, "chr2_74Mb_wideLots"],
              ["chr2" ,157040000, 157220000, "chr2_157Mb"],
              ["chr2" ,160622000, 161132000, "chr2_161Mb"],
              ["chr2" ,160712000, 160952000, "chr2_160Mb"],
              ["chr4" , 33203000,  33253000, "chr4_33Mb"],
              ["chr4" , 59537000,  59637000, "chr4_59Mb"],
              ["chr4" ,103104000, 103224000, "chr4_103Mb"],
              ["chr4" ,123655000, 123755000, "chr4_123Mb"],
              ["chr4" ,126961000, 127031000, "chr4_127Mb"],
              ["chr6" , 52135120,  52282912, "chr6_52Mb_HoxA"],
              ["chr6" , 53000000,  53350000, "chr6_53Mb"],
              ["chr6" , 72350449,  72361990, "chr6_72Mb"],
              ["chr7" ,100842558 ,100860108, "chr7_100Mb_wideNoOri"],
              ["chr8" , 70461513,  70618923, "chr8_70Mb"],
              ["chr8" , 79251596,  79688107, "chr8_79Mb"],
              ["chr8" , 85382908,  85746511, "chr8_85Mb"],
              ["chr8" ,106023193, 106246403, "chr8_106Mb"],
              ["chr8" ,106854116, 107077326, "chr8_107Mb"],
              ["chr8" ,109681099, 109904309, "chr8_109Mb"],
              ["chr8" ,122452198, 122572643, "chr8_122Mb"],
              ["chr11", 75700000,  75825000, "chr11_75Mb"],
              ["chr11", 77439474,  77540385, "chr11_77Mb_wideLots"],
              ["chr15", 74487685,  74629350, "chr15_74Mb_complex"],
              ["chr15", 84668915,  84707234, "chr15_84Mb_1_0_twoWide"],
              ["chr18", 61672481,  61700287, "chr18_61Mb_wideNoOri"])
      .into {slices2use; s2u}

process makeSlices {

  tag {sName}

  publishDir params.outdirSlices, mode: 'copy', overwrite: true, pattern: '*.tab'
  publishDir params.outdirSlices, mode: 'copy', overwrite: true, pattern: '*.bed'
  publishDir params.outdirSlices, mode: 'copy', overwrite: true, pattern: '*.bedgraph'

  input:
  file (af) from af10.collect()
  set val(nCS), val(nFrom), val(nTo), val(sName) from slices2use
  file(bed)                                      from oriSSDSBEDs_mm10_4Slices.collect()
  file(win4bg)                                   from winFiles_a.collect()
  file(ori)                                      from allOriginsBED_b.collect()
  file(tss)                                      from refseqGeneBED
  file(gtf)                                      from gencodeGTFa
  file(cpg)                                      from mm10CpGpeaksBED_a
  file(cpgDensityBG)                             from mm10CpGdensityBG_a

  output:
  //file("chr*DZ.tab")       into sliceDZ
  file("*ss*bedgraph")         into (sliceBG,sliceBG_a,sliceBG_b)
  file("chr*all.dz.tab")       into (sliceTab,sliceTab_a,sliceTab_b)
  file("chr*ori.bed")          into (sliceOri,sliceOri_a,sliceOri_b)
  file("chr*tss.bed")          into (sliceTSS,sliceTSS_a,sliceTSS_b)
  file("chr*ATACSeq.bed")      into (sliceATAC,sliceATAC_a,sliceATAC_b)
  file("chr*CGI.bed")          into (sliceCGI,sliceCGI_a,sliceCGI_b)
  file("chr*genes.bed")        into (sliceGenes,sliceGenes_a,sliceGenes_b)
  file("chr*g4.bed")           into (sliceG4,sliceG4_a,sliceG4_b)
  file("chr*CpGpeaks.bed")     into (sliceCpGpeaks,sliceCpGpeaks_a,sliceCpGpeaks_b)
  file("chr*CpGdens.bedgraph") into (sliceCpGdensity,sliceCpGdensity_a,sliceCpGdensity_b)

  script:
  """
  echo -e "${nCS}\\t${nFrom}\\t${nTo}"  >region.bed

  ori="hiconf_origins.mm10.bedgraph"

  intersectBed -a ${tss} -b region.bed -wa -u |cut -f1-6 >${sName}.tss.bed
  intersectBed -a \$ori  -b region.bed -wa -u |cut -f1-3 |sort -k1,1 -k2n,2n >${sName}.ori.bed

  intersectBed -a ${params.datadir}/g4/mm10.g4.bed -b region.bed >g4.bed
  perl -lane \'@X=split(":",\$F[4]); \$from=(\$F[5] eq "+")?\$F[1]:\$F[2]; \$to=(\$F[5] eq "+")?\$F[2]:\$F[1]; print join("\\t",\$F[0],\$from,\$to,@F[3..5]) if (\$X[0] >= 6)\' g4.bed |cut -f1-6 >${sName}.g4.bed

  intersectBed -a ${gtf} -b region.bed |grep -P \'\\stranscript\\s\' |cut -f1,4,5,7 |perl -lane \'print join("\\t",@F[0..2],\'x\',\'x\',\$F[3])\' |sort -k1,1 -k2n,2n >genes.bed
  mergeBed -i genes.bed -s -c 4,5,6 -o distinct,distinct,distinct >${sName}.genes.bed

  intersectBed -a ${params.annotationdir}/CGI.mm10.bedgraph     -b region.bed |cut -f1-3 |sort -k1,1 -k2n,2n -k3n,3n >${sName}.CGI.bed
  intersectBed -a ${params.annotationdir}/ATACSeq.mm10.bedgraph -b region.bed |cut -f1-3 |sort -k1,1 -k2n,2n -k3n,3n >${sName}.ATACSeq.bed

  intersectBed -a ${cpg}          -b region.bed -wa -u >${sName}.CpGpeaks.bed
  intersectBed -a ${cpgDensityBG} -b region.bed -wa -u >${sName}.CpGdens.bedgraph

  intersectBed -a win2ks100.bed -b region.bed -wa -u >wins.${sName}.w2000s100.bed
  intersectBed -a win1ks147.bed -b region.bed -wa -u >wins.${sName}.w1000s147.bed
  intersectBed -a win500s50.bed -b region.bed -wa -u >wins.${sName}.w500s50.bed
  intersectBed -a win500.bed    -b region.bed -wa -u >wins.${sName}.ws500.bed

  for ssBED in *type1.bed; do

    outBG2k100f=\${ssBED/bed/${sName}.F.w2000s100.bedgraph}
    outBG2k100r=\${ssBED/bed/${sName}.R.w2000s100.bedgraph}
    outBG2k100t=\${ssBED/bed/${sName}.T.w2000s100.bedgraph}

    outBG1k147f=\${ssBED/bed/${sName}.F.w1000s147.bedgraph}
    outBG1k147r=\${ssBED/bed/${sName}.R.w1000s147.bedgraph}
    outBG1k147t=\${ssBED/bed/${sName}.T.w1000s147.bedgraph}

    outBG500s50f=\${ssBED/bed/${sName}.F.w500s50.bedgraph}
    outBG500s50r=\${ssBED/bed/${sName}.R.w500s50.bedgraph}
    outBG500s50t=\${ssBED/bed/${sName}.T.w500s50.bedgraph}

    outBG500f=\${ssBED/bed/${sName}.F.ws500.bedgraph}
    outBG500r=\${ssBED/bed/${sName}.R.ws500.bedgraph}
    outBG500t=\${ssBED/bed/${sName}.T.ws500.bedgraph}

    nSS=`cat \$ssBED |wc -l`

    intersectBed -a \$ssBED -b region.bed |sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 |uniq >ss.bed

    grep -P  \'\\+\' ss.bed >FWD.bed
    grep -vP \'\\+\' ss.bed >REV.bed

    intersectBed -a wins.${sName}.w2000s100.bed -b FWD.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG2k100f
    intersectBed -a wins.${sName}.w2000s100.bed -b REV.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG2k100r
    intersectBed -a wins.${sName}.w2000s100.bed -b  ss.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG2k100t

    intersectBed -a wins.${sName}.w1000s147.bed -b FWD.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG1k147f
    intersectBed -a wins.${sName}.w1000s147.bed -b REV.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG1k147r
    intersectBed -a wins.${sName}.w1000s147.bed -b  ss.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG1k147t

    intersectBed -a wins.${sName}.w500s50.bed   -b FWD.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500s50f
    intersectBed -a wins.${sName}.w500s50.bed   -b REV.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500s50r
    intersectBed -a wins.${sName}.w500s50.bed   -b  ss.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500s50t

    intersectBed -a wins.${sName}.ws500.bed     -b FWD.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500f
    intersectBed -a wins.${sName}.ws500.bed     -b REV.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500r
    intersectBed -a wins.${sName}.ws500.bed     -b  ss.bed -c |awk '{print \$1, \$2, \$3, 1000000*\$4/('\$nSS')}' |sort -k1,1 -k2n,2n >\$outBG500t

  done

  echo -e \'cs\\tpos\\tcover\\tname\\tws\\tfr\' >${sName}.all.dz.tab

  for ssBG in *bedgraph; do
    dz=\${ssBG/bedgraph/DZ.tab}
    nm=\${ssBG/bedgraph/DZ.tab}
    perl -lane \'\$min = int((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$min,\$F[3])\' \$ssBG >\$dz
    perl -lane \'\$bg="\'\$ssBG\'"; \$bg =~ /^(.+)\\.oriSSDS.+\\.([RTF])\\.(\\S+?)\\./; (\$nm,\$FR, \$ws) = (\$1,\$2,\$3); \$min = int((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$min,\$F[3],\$nm,\$ws,\$FR)\' \$ssBG >>${sName}.all.dz.tab
  done
  """
  }

process drawRNAhydrolysisFig {

  publishDir params.outdirFigs, mode: 'copy', overwrite: true

  input:
  file (af) from af11.collect()
  file(sliceTab)   from sliceTab_b.collect()
  file(sliceOri)   from sliceOri_b.collect()
  file(sliceTSS)   from sliceTSS_b.collect()
  file(sliceATAC)  from sliceATAC_b.collect()
  file(sliceCGI)   from sliceCGI_b.collect()
  file(sliceGenes) from sliceGenes_b.collect()
  file(sliceG4)    from sliceG4_b.collect()
  file(oriTable)   from allOriginRData_a.collect()
  file(dtMatrix)  from oriSSDS_dtMatrix_b.collect()

  output:
  file('Pratto*png')   into suppFigRNAhPNG
  file('Pratto*svg')   into suppFigRNAhSVG
  file('Pratto*pdf')   into suppFigRNAhPDF

  script:
  """
  for z in *gz; do
    unzipped=\${z/.gz/}
    gunzip -c \$z > \$unzipped
  done

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  R --no-save <accessoryFiles/scripts/R/plotRNAhydrolysisFig.R ||true
  """
  }

process drawSNSPeakCallingProblemsFig {

  publishDir params.outdirFigs, mode: 'copy', overwrite: true

  input:
  file (af) from af12.collect()

  output:
  file('Pratto*png')   into suppFigSNSpeaksPNG
  file('Pratto*svg')   into suppFigSNSpeaksSVG
  file('Pratto*pdf')   into suppFigSNSpeaksPDF

  script:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134988/suppl/GSE134988%5FsiNC%2DNS%5Fpeaks%2Ebed%2Egz
  gunzip GSE134988_siNC-NS_peaks.bed.gz

  slopBed -i GSE134988_siNC-NS_peaks.bed -l -0.5 -r -0.5 -pct -g ${hg38IDX} |grep -P '^chr' |slopBed -i - -l 1500 -r 1500 -g ${hg38IDX} |sort -k1,1 -k2n,2n >long.bed
  mergeBed -i long.bed -c 3 -o count |perl -lane 'print join("\\t","Long_2020",\$F[3]) if (\$F[2]-\$F[1] < 10000)' >snsClusters.tab

  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68347/suppl/GSE68347%5FInitiation%5FSites%2EbedGraph%2Egz
  gunzip GSE68347_Initiation_Sites.bedGraph.gz
  cut -f1-3 GSE68347_Initiation_Sites.bedGraph >GSE68347_Initiation_Sites.mm9.bed

  wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz
  liftOver GSE68347_Initiation_Sites.mm9.bed  mm9ToMm10.over.chain.gz  GSE68347_Initiation_Sites.mm10.bed na

  slopBed -i GSE68347_Initiation_Sites.mm10.bed -l -0.5 -r -0.5 -pct -g ${mm10IDX} |slopBed -i - -l 1500 -r 1500 -g ${mm10IDX} |sort -k1,1 -k2n,2n >cayrou.bed
  mergeBed -i cayrou.bed -c 3 -o count |perl -lane 'print join("\\t","Cayrou_2015",\$F[3]) if (\$F[2]-\$F[1] < 10000)' >>snsClusters.tab

  ln -s accessoryFiles/img/OriSSDS_v_SNSSeq_Mechali_oneOrigin.png .

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  R --no-save <accessoryFiles/scripts/R/plotSNSpeakClusteringFig.R ||true
  """
  }


process getCpGvOriPredictiveAbility {

  publishDir params.outdirGomez, mode: 'copy', overwrite: true

  input:
  file (af) from af28.collect()
  file(winz) from winFiles_f.collect()
  file(ori) from mouseOriginsBG_h

  output:
  file 'CpG_overlapStats.txt'      into mm10CpGOLStats

  script:
  """
  ## Get CpG density in 1 Kb wins 100bp slide
  perl ${params.codedir}/getDiNTsFromFA.pl ${mm10FA} CG >mm10_CpG.bed

  mapBed -a win1ks100.bed -b mm10_CpG.bed -g ${mm10IDX} -c 1 -o count |perl -lane 'print join("\\t",\$F[0],\$F[1]+450,\$F[2]-450,\$F[3])' >mm10_CpG.w1ks100.tmp
  sort -k1,1 -k2n,2n -k3n,3n mm10_CpG.w1ks100.tmp >mm10_CpG.w1ks100.bedgraph

  ## Get overlap stats of CpG peaks with Oris
  echo -e "pcGC\\tatOri\\ttotGC\\ttotOri" >CpG_overlapStats.txt

  cp mm10_CpG.w1ks100.bedgraph mm10_CpG.w1ks100.forT.bedgraph

  for p in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15; do

    perl -lane '\$cpg = \$F[3]/1000; print \$_ if (\$cpg > '\$p')' mm10_CpG.w1ks100.forT.bedgraph >CpG_thresholded.bed

    ## Speed this up ...
    cp CpG_thresholded.bed mm10_CpG.w1ks100.forT.bedgraph

  	slopBed   -i CpG_thresholded.bed           -g ${mm10IDX} -b 100               |sort -k1,1 -k2n,2n              >CpG_thresholded.slop100.bed

    mergeBed  -i CpG_thresholded.slop100.bed                                      |sort -k1,1 -k2n,2n              >CpG_thresholded.slopMerge.bed

    slopBed   -i CpG_thresholded.slopMerge.bed -g ${mm10IDX} -pct -r -0.5 -l -0.5 |slopBed -i - -g accessoryFiles/genomeFiles//mm10_genome.fa.fai -b 150  |sort -k1,1 -k2n,2n              >mm10_CpG_peaks.\$p.bed

    ol=`intersectBed -a ${ori}  -b "mm10_CpG_peaks."\$p".bed" -u |grep -vP 'chr[XYM]' |wc -l`

    tot=`cat "mm10_CpG_peaks."\$p".bed" |wc -l`

    ori=`cat ${ori}|wc -l`

    echo -e "\$p\\t\$ol\\t\$tot\\t\$ori" >>CpG_overlapStats.txt

  done
  """
  }

process getTFBSfiles {

  input:
  file (hs) from hotspotBG500
  file (af) from af41.collect()

  output:
  path('TFBS_*') into TFBS
  file('*jar')   into apps
  file('HOCO*')  into motifs
  path('pwm')    into pwm
  path('pvals')  into pvals

  script:
  """
  wget http://opera.autosome.ru/downloads/ape-3.0.2.jar
  wget https://raw.githubusercontent.com/autosome-ru/sarus/master/releases/sarus-2.0.1.jar
  wget http://gtrd.biouml.org/downloads/19.10/chip-seq/Mus%20musculus_meta_clusters.interval.gz -O TF_peaks.tar.gz
  wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/MOUSE/mono/HOCOMOCOv11_core_pwm_MOUSE_mono.tar.gz

  #cp accessoryFiles/TF/*jar .
  #cp accessoryFiles/TF/*gz .

  zcat TF_peaks.tar.gz | grep -v uniprot | grep -P 'chr[0-9]+\\s' |awk -F"\\t" '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"toupper(\$6)"\\t+">toupper(\$6)".TFBS.bed"}'

  ## Use DSB hotspots +- 250 bp for PRDM9
  grep -P 'chr[0-9]+\\s' ${hs} |perl -lane 'print join("\\t",@F[0..3],"PRDM9","+")' >PRDM9.TFBS.bed

  for tfbs in *.TFBS.bed; do
    sort -k1,1 -k2n,2n \$tfbs -o \$tfbs
  done

  tar -zxvf HOCOMOCOv11_core_pwm_MOUSE_mono.tar.gz
  java -cp ape-3.0.2.jar ru.autosome.ape.PrecalculateThresholds pwm pvals --background .3,.2,.2,.3

  ## Rename PWMs
  for p in pwm/*pwm; do
    t1=\${p/pwm\\//};
    tf=\${t1/_MOUSE*pwm/};
    if [ -e \$tf".TFBS.bed" ]; then
      mv \$p pwm/\$tf".pwm"
    else
      rm \$p
    fi
  done

  ## Rename Threshold Files
  for t in pvals/*thr; do
    t1=\${t/pvals\\//}
    tf=\${t1/_MOUSE*thr/}
    if [ -e \$tf".TFBS.bed" ]; then
      mv \$t pvals/\$tf".thr"
    else
      rm \$t
    fi
  done

  for t in *.TFBS.bed; do
    tf=\${t/.TFBS.bed/}
    if [ ! -e 'pwm/'\$tf'.pwm' ]; then
      rm \$t
    fi
  done

  ## Split to subfolders
  i=1
  while read l; do mkdir -p TFBS_\$i;cp \$l TFBS_\$((i++));done< <(ls *.TFBS.bed|xargs -n10)
  """
  }

process makeROCForTF {

  publishDir "${params.outdirFigs}/TF_ROC", mode: 'copy', overwrite: true

  input:
  file (af) from af42.collect()
  each path(TFBS)
  file(apps)
  file(motifs)
  path(pwm)
  path(pvals)

  output:
  file('*ROC.tab') into ROCdata
  file('*png')     into png
  file('*pdf')     into pdf

  script:
  """

  for f in ${TFBS}/* ; do
    ln -s \$f .
  done

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  for tfbed in *TFBS.bed; do
    tf=\${tfbed/.TFBS.bed/}
    tfbed=\$tf".TFBS.bed"
    pwm="pwm/"\$tf".pwm"
    mbed=\$tf".motifs.bed"
    pktab=\$tf".peaks_with_motifs.bed"

    score_thresh=`perl -lane 'print \$F[0]; exit if (\$F[1] < 0.005)' pvals/\$tf".thr" |tail -n1`

    npeaks=`cat \$tfbed |wc -l`

    if [[ \$npeaks -gt 12000 ]]; then
    	java -jar sarus-2.0.1.jar ${mm10FA} \$pwm \$score_thresh --output-bed | intersectBed -a - -b \$tfbed -c >\$mbed
    	sort -k1,1 -k2n,2n \$mbed -o \$mbed
    	intersectBed -a \$tfbed -b \$mbed -wao -sorted |slopBed -l -0.5 -r -0.5 -pct -g ${mm10IDX} -i - |perl -pi -e 's/\\t\\./\\t99999/g' |mergeBed -i - -c 11 -o min |perl -lane '\$F[3] =~ s/99999/-999/; print join("\\t",@F)' >\$pktab

    	ln \$mbed motifs.bed
    	ln \$pktab peaks_with_motifs.bed

    	R --no-save <accessoryFiles/scripts/R/ROC_from_TFBS.R

    	unlink motifs.bed
    	unlink peaks_with_motifs.bed
    else
     	echo "SKIP \$tfbed"
    fi
  done
  """
  }

process originsVCpGs {

  publishDir params.outdirImages,  mode: 'copy', overwrite: true, pattern: '*p??'
  publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*tab'

  input:
  file (af)      from af13.collect()
  file(ROC)      from ROCdata.collect()
  file(ori)      from mouseOriginsBG_b
  file(tss)      from tssNoMergeBED1bp
  file(oriRdata) from mouseOriginsHCRdata
  file(mm10CpGOLStats)

  output:
  file("ori_v_cgis.tab")               into oriVCpGsTAB
  file("*png")                         into oriVCpGsPNG
  file("*pdf")           optional true into oriVCpGsPDF

  script:
  """
  echo 'test'
  mapBed -a ${ori} -b ${params.annotationdir}/CGI.mm10.bedgraph            -c 4 -o sum >oriCGI.tab
  mapBed -a ${ori} -b ${params.annotationdir}/CGI_at_ATACSeq.mm10.bedgraph -c 4 -o sum >oriCGIatATAC.tab
  mapBed -a ${ori} -b ${params.annotationdir}/CGI_no_ATACSeq.mm10.bedgraph -c 4 -o sum >oriCGInoATAC.tab

  echo -e "typeA\\ttypeB\\tATot\\tBTot\\tA\\tB\\tAB\\tBA\\tpcAB\\tpcBA" >ori_v_cgis.tab

  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/CGI_mm10.bedgraph            --nA ori --nB allCGIs        >>ori_v_cgis.tab
  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/CGI_at_ATACSeq.mm10.bedgraph --nA ori --nB CGI_at_ATACSeq >>ori_v_cgis.tab
  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/CGI_no_ATACSeq.mm10.bedgraph --nA ori --nB CGI_no_ATACSeq >>ori_v_cgis.tab
  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/ATACSeq.mm10.bedgraph        --nA ori --nB allATACSeq     >>ori_v_cgis.tab
  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/ATACSeq_at_CGI.mm10.bedgraph --nA ori --nB ATACSeq_at_CGI >>ori_v_cgis.tab
  perl ${params.codedir}/beds2Venn.pl --a ${ori} --b ${params.annotationdir}/ATACSeq_no_CGI.mm10.bedgraph --nA ori --nB ATACSeq_no_CGI >>ori_v_cgis.tab

  cut -f1-3 ${ori} >ori.bed
  intersectBed -a ${ori}                                        -b ${params.annotationdir}/CGI_mm10.bedgraph     -c >x.x
  intersectBed -a ${ori}                                        -b ${params.annotationdir}/ATACSeq.mm10.bedgraph -c >y.y
  paste x.x y.y |perl -lane 'print join("\\t","ori",(\$F[4]>0?"TRUE":"FALSE"),(\$F[9]>0?"TRUE":"FALSE"))' >ori.ol.tab

  intersectBed -a ${params.annotationdir}/ATACSeq.mm10.bedgraph -b ${ori}                                        -c >a.x
  intersectBed -a ${params.annotationdir}/ATACSeq.mm10.bedgraph -b ${params.annotationdir}/CGI_mm10.bedgraph     -c >a.y
  paste a.x a.y |perl -lane 'print join("\\t","ATAC",(\$F[4]>0?"TRUE":"FALSE"),(\$F[9]>0?"TRUE":"FALSE"))' >ATAC.ol.tab

  intersectBed -a ${params.annotationdir}/CGI_mm10.bedgraph     -b ${ori}                                        -c >c.x
  intersectBed -a ${params.annotationdir}/CGI_mm10.bedgraph     -b ${params.annotationdir}/ATACSeq.mm10.bedgraph -c >c.y
  paste c.x c.y |perl -lane 'print join("\\t","CGI",(\$F[4]>0?"TRUE":"FALSE"),(\$F[9]>0?"TRUE":"FALSE"))' >CGI.ol.tab

  cat *ol.tab >all.ol.tab

  intersectBed -a ${params.annotationdir}/CGI_mm10.bedgraph -b ${params.annotationdir}/ATACSeq.mm10.bedgraph -wo |perl -lane 'print join("\t",\$F[0],\$F[1]<\$F[5]?\$F[1]:\$F[5],\$F[2]>\$F[6]?\$F[2]:\$F[6])' |uniq >ATAC-CGI.bed
  intersectBed -a ${params.annotationdir}/CGI_mm10.bedgraph -b ${params.annotationdir}/ATACSeq.mm10.bedgraph -v |cut -f1-3 >CGI_no_ATACSeq.bed
  intersectBed -b ${params.annotationdir}/CGI_mm10.bedgraph -a ${params.annotationdir}/ATACSeq.mm10.bedgraph -v |cut -f1-3 >ATACSeq_no_CGI.bed

  olATACCGI=`intersectBed  -a ATAC-CGI.bed                                  -b hiconf_origins.mm10.bedgraph -u |wc -l`
  olCGI=`intersectBed      -a CGI_no_ATACSeq.bed                            -b hiconf_origins.mm10.bedgraph -u |wc -l`
  olATAC=`intersectBed     -a ATACSeq_no_CGI.bed                            -b hiconf_origins.mm10.bedgraph -u |wc -l`
  olCGIall=`intersectBed   -a ${params.annotationdir}/CGI_mm10.bedgraph     -b hiconf_origins.mm10.bedgraph -u |wc -l`
  olATACall=`intersectBed  -a ${params.annotationdir}/ATACSeq.mm10.bedgraph -b hiconf_origins.mm10.bedgraph -u |wc -l`

  totATACCGI=`cat  ATAC-CGI.bed                                 |wc -l`
  totCGI=`cat CGI_no_ATACSeq.bed                                |wc -l`
  totATAC=`cat ATACSeq_no_CGI.bed                               |wc -l`
  totCGIall=`cat ${params.annotationdir}/CGI_mm10.bedgraph      |wc -l`
  totATACall=`cat ${params.annotationdir}/ATACSeq.mm10.bedgraph |wc -l`

  echo -e 'type\tatOri\ttot'                      >olStat.tab
  echo -e "CGI\t\$olCGI\t\$totCGI"                >>olStat.tab
  #echo -e "CGIall\t\$olCGIall\t\$totCGIall"      >>olStat.tab
  echo -e "ATAC-Seq peaks\t\$olATAC\t\$totATAC"   >>olStat.tab
  #echo -e "ATACall\t\$olATACall\t\$totATACall"   >>olStat.tab
  echo -e "ATAC & CGI\t\$olATACCGI\t\$totATACCGI" >>olStat.tab

  ## Add globals to .Renviron (so we can see them inside R!)

  ## Next plots
  bedtools slop -i ${params.annotationdir}/ATACSeq.mm10.bedgraph -r -0.5 -l -0.5 -pct -g ${mm10IDX} |cut -f1-3 |sort -k1,1 -k2n,2n >atacPP.bed
  bedtools slop -i ${params.annotationdir}/CGI_mm10.bedgraph     -r -0.5 -l -0.5 -pct -g ${mm10IDX} |cut -f1-3 |sort -k1,1 -k2n,2n >cgiPP.bed
  bedtools slop -i ${tss}                                        -r -0.5 -l -0.5 -pct -g ${mm10IDX} |cut -f1-3 |sort -k1,1 -k2n,2n >tssPP.bed
  bedtools slop -i hiconf_origins.mm10.bedgraph                  -r -0.5 -l -0.5 -pct -g ${mm10IDX} |cut -f1-3 |sort -k1,1 -k2n,2n >oriPP.bed

  bedtools shuffle -g ${mm10IDX} -seed 42 -i oriPP.bed  |sort -k1,1 -k2n,2n >oriRPP.bed

  rm -f closest.tab
  bedtools closest -a oriPP.bed  -b atacPP.bed -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Origins","atac")' >>closest.tab
  bedtools closest -a oriPP.bed  -b cgiPP.bed  -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Origins","cgi")'  >>closest.tab
  bedtools closest -a oriPP.bed  -b tssPP.bed  -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Origins","tss")'  >>closest.tab
  bedtools closest -a oriRPP.bed -b atacPP.bed -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Random" ,"atac")' >>closest.tab
  bedtools closest -a oriRPP.bed -b cgiPP.bed  -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Random" ,"cgi")'  >>closest.tab
  bedtools closest -a oriRPP.bed -b tssPP.bed  -d -t first  |cut -f 7 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Random" ,"tss")'  >>closest.tab

  perl -lane 'print join("\t",@F,join("_",@F))' oriPP.bed  >oriPP.bed4
  perl -lane 'print join("\t",@F,join("_",@F))' oriRPP.bed >oriRPP.bed4

  bedtools closest -N -a oriPP.bed4  -b oriPP.bed4  -d -t first |cut -f 9 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Origins","Origins")' >oriVori.tab
  bedtools closest -N -a oriRPP.bed4 -b oriRPP.bed4 -d -t first |cut -f 9 |perl -lane 'print join("\\t",sprintf("%i",\$_),"Random","Random")'  >>oriVori.tab

  # Final ones
  bedtools slop -l 5500 -r 5500 -g ${mm10IDX} -i oriPP.bed |sort -k1,1 -k2n,2n >ori11k.bed

  intersectBed -a ori11k.bed -b accessoryFiles/data/g4/mm10.qparser2.g4.bed -wao |\
  perl -lane 'next if (\$_ =~ /\\s\\-1\\s+\\-1/); \$mid=((\$F[2]+\$F[1])/2); \$g4=((\$F[4]+\$F[5])/2); \$d=\$g4-\$mid; print join("\\t",\$mid,\$g4,\$d,\$F[8])' >g4_vs_origins.tab

  ########################################################################################
  cp ${params.codedir}/CalcPhysicoChemicalPropsFromFA.pl .

  bedtools slop -l -0.5 -r -0.5 -pct -i ${ori} -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |cut -f1-3 >oriMMnt.bed

  bedtools slop -l 10000 -r 10000 -i oriMMnt.bed -g ${mm10IDX} >oriMM.10k.bed

  bedtools getfasta -fi ${mm10FA} -bed oriMM.10k.bed -fo oriMM.10k.fa

  perl CalcPhysicoChemicalPropsFromFA.pl oriMM.10k.fa 25 >mm10Origins.bendability.txt
  ########################################################################################

  # Gather TFBS ROC data #################################################################
  cat *ROC.tab |grep sens |head -n1  >ROCdata.tab
  cat *ROC.tab |grep -v sens        >>ROCdata.tab
  ########################################################################################

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  R --no-save <accessoryFiles/scripts/R/plot_mm10_Ori_v_CGIs.R ||true

  """
  }

process analyzeOriClusters {

  publishDir params.outdirImages,  mode: 'copy', overwrite: true, pattern: '*p??'

  input:
  file (af) from af14.collect()
  file(ori)      from mouseOriginsBG_e
  file(bw)       from oriSSDSBW1.collect()
  file(winz)     from winFiles_b.collect()
  file(bgCpG)    from mm10CpGdensityBG_b
  file(cpgPeaks) from mm10CpGpeaksBED_b
  file(gn)       from sliceGenes_a.collect()
  file(at)       from sliceATAC_a.collect()
  file(ta)       from sliceTab_a.collect()
  file(or)       from sliceOri_a.collect()
  file(ts)       from sliceTSS_a.collect()
  file(cg)       from sliceCGI_a.collect()
  file(g4)       from sliceG4_a.collect()
  file(sliceBG)  from sliceBG_a.collect()
  file(cpgPeak)  from sliceCpGpeaks_a.collect()
  file(cpgDens)  from sliceCpGdensity_a.collect()

  output:
  file("*png")           into oriCustersPNG
  file("*pdf")           into oriClustersPDF

  script:
  """
  cut -f1-3 ${ori} |sort -k1,1 -k2n,2n -k3n,3n >origins.bed
  perl -lane 'print \$_ if ((\$F[2]-\$F[1]) >= 10000)' origins.bed |sort -k1,1 -k2n,2n >ori10k.bed
  perl -lane 'print \$_ if ((\$F[2]-\$F[1]) <= 7000)'  origins.bed |sort -k1,1 -k2n,2n >oriNarrow.bed

  ## All origins first
  #intersectBed -a origins.bed -b mm10_CpG.w1ks100.bedgraph -wao -sorted |perl -lane '\$cpg = \$F[6]/1000; print \$_ if (\$cpg > 0.05)' |cut -f4-7 >CpG_atOriAll.bed
  #slopBed   -i CpG_atOriAll.bed           -g ${mm10IDX} -b 100                                                   |sort -k1,1 -k2n,2n              >CpG_atOriAll.slop100.bed
  #mergeBed  -i CpG_atOriAll.slop100.bed                                                                          |sort -k1,1 -k2n,2n              >CpG_atOriAll.slopMerge.bed
  #slopBed   -i CpG_atOriAll.slopMerge.bed -g ${mm10IDX} -pct -r -0.5 -l -0.5 |slopBed -i - -g ${mm10IDX} -b 150  |sort -k1,1 -k2n,2n              >CpG_atOriAll.peaks.bed

  intersectBed -a mm10_CpG_peaks.bed -b origins.bed -u >CpG_atOriAll.peaks.bed

  computeMatrix reference-point --referencePoint center -a 5000 -b 5000 \
    -R CpG_atOriAll.peaks.bed \
    -S mm10_WT_Rep2.NEG.frags.RPKM.bigwig mm10_WT_Rep2.POS.frags.RPKM.bigwig \
      mm10_sperm_Rep1.NEG.frags.RPKM.bigwig mm10_sperm_Rep1.POS.frags.RPKM.bigwig \
    -p max -o CpG_atOriAll.peaks.Matrix.gz \
    --binSize 100
    #mm10_WT_Rep1.NEG.frags.RPKM.bigwig mm10_WT_Rep1.POS.frags.RPKM.bigwig \
    #mm10_WT_Rep3.NEG.frags.RPKM.bigwig mm10_WT_Rep3.POS.frags.RPKM.bigwig \

  plotHeatmap --matrixFile CpG_atOriAll.peaks.Matrix.gz \
    -o originZones.CpG_atOriAll.peaks.deeptools.png \
    --xAxisLabel "Distance to CpG peak (Kb)" \
    --colorMap Reds \
    --samplesLabel "Crick (Testis)" "Watson (Testis)" "Crick (Sperm)" "Watson (Sperm)" \
    --dpi 300 \
    --regionsLabel "CpG peaks" \
    --refPointLabel CpG --averageTypeSummaryPlot mean
    #--samplesLabel "Crick (Rep1)" "Watson (Rep1)" "Crick (Rep2)" "Watson (Rep2)" "Crick (Rep3)" "Watson (Rep3)" "Crick (Sperm)" "Watson (Sperm)" \

  ## Now 10k origins
  intersectBed -a mm10_CpG_peaks.bed -b ori10k.bed -u >CpG_atOri10k.peaks.bed

  computeMatrix reference-point --referencePoint center -a 5000 -b 5000 \
    -R CpG_atOri10k.peaks.bed \
    -S mm10_WT_Rep2.NEG.frags.RPKM.bigwig mm10_WT_Rep2.POS.frags.RPKM.bigwig \
      mm10_sperm_Rep1.NEG.frags.RPKM.bigwig mm10_sperm_Rep1.POS.frags.RPKM.bigwig \
    -p max -o CpG_atOri10k.peaks.Matrix.gz \
    --binSize 100
    #mm10_WT_Rep1.NEG.frags.RPKM.bigwig mm10_WT_Rep1.POS.frags.RPKM.bigwig \
    #mm10_WT_Rep3.NEG.frags.RPKM.bigwig mm10_WT_Rep3.POS.frags.RPKM.bigwig \

  plotHeatmap --matrixFile CpG_atOri10k.peaks.Matrix.gz\
    -o originZones.CpG_atOri10k.peaks.deeptools.png \
    --xAxisLabel "Distance to CpG peak (Kb)" \
    --colorMap Reds \
    --samplesLabel "Crick (Testis)" "Watson (Testis)" "Crick (Sperm)" "Watson (Sperm)" \
    --dpi 300 \
    --regionsLabel "CpG peaks" \
    --refPointLabel CpG --averageTypeSummaryPlot mean
  #--samplesLabel "Crick (Rep1)" "Watson (Rep1)" "Crick (Rep2)" "Watson (Rep2)" "Crick (Rep3)" "Watson (Rep3)" "Crick (Sperm)" "Watson (Sperm)" \

  gunzip CpG_atOriAll.peaks.Matrix.gz
  gunzip CpG_atOri10k.peaks.Matrix.gz

  ## Get CpG distrib v origin center
  slopBed -i oriNarrow.bed -pct -b -0.5 -g ${mm10IDX} |slopBed -i - -b 3000 -g ${mm10IDX} |sort -k1,1 -k2n,2n >ori6k.bed
  intersectBed -a ori6k.bed  -b CpG_atOriAll.peaks.bed -c   |grep -P '\\s1\$' |cut -f1-3 >oneCpg.bed
  intersectBed -a oneCpg.bed -b CpG_atOriAll.peaks.bed -wao |perl -lane 'use Math::Round; \$mo=round((\$F[1]+\$F[2])/2); \$mc=round((\$F[4]+\$F[5])/2); print (\$mo-\$mc)' >oriVCpG.tab

  # Count CpG peaks per Zone
  intersectBed -a origins.bed -b CpG_atOriAll.peaks.bed -c >CpGpeaksAtOrigins.txt

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  R --no-save <accessoryFiles/scripts/R/drawZonesFigure.R ||true

  """
  }

process makeFigure1 {

  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*igu*.png'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*igu*.pdf'

  input:
  file (af) from af15.collect()
  file(oriTable) from allOriginRData.collect()
  file(dtMatrix) from oriSSDS_dtMatrix_a.collect()
  file(gn)       from sliceGenes.collect()
  file(at)       from sliceATAC.collect()
  file(ta)       from sliceTab.collect()
  file(or)       from sliceOri.collect()
  file(ts)       from sliceTSS.collect()
  file(cg)       from sliceCGI.collect()
  file(g4)       from sliceG4.collect()
  file(sliceBG)  from sliceBG.collect()

  output:
  file("*igure1*png")             into figure1PNG
  file("*igure1*pdf")             into figure1PDF

  script:
  """
  cp accessoryFiles/img/oriSSDSHSlim.png .

  gunzip mm10_WT_Rep2.origins.fragsMatrix.gz -c >mm10_WT_Rep2.origins.fragsMatrix

  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow    ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"     >>./.Renviron

  R --no-save <accessoryFiles/scripts/R/drawFigure1.R
  """
  }

process makeMouseRTSeqTable {

  tag {sName}

  publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*tab'

  input:
  file (af) from af16.collect()

  output:
  file("RT.mm10.tab")               into mm10RTtable

  script:
  """
  for f in accessoryFiles/RTSeq/bedgraph/final/raw/*.mm10.RT.bedgraph; do
    n=`basename \$f`
    nm=\${n/.mm10.RT.bedgraph/}
    echo \$nm >\$nm.RT
    cut -f4 \$f >>\$nm.RT
    echo -e "cs\\tfrom\\tto" >pos.bed
    cut -f1-3 \$f >>pos.bed
  done

  paste pos.bed *RT >RT.mm10.tab
  """
  }

// if (params.getStatsAgain){
//   process getModelGridSearchResults {
//     scratch '/lscratch/$SLURM_JOBID'
//     clusterOptions ' --gres=lscratch:40'
//     echo true
//     cpus 1
//     memory '128g'
//
//     time { 1.hour }
//
//     module 'bedtools/2.27.1'
//     module 'ucsc/388'
//     module 'R/3.5.2'
//
//     publishDir params.outdirRTables, mode: 'copy', overwrite: true
//
//     input:
//
//     output:
//     file("stats.mm10.tab")               into mouseModellingStats
//
//     script:
//     """
//     ## Add globals to .Renviron (so we can see them inside R!)
//
//     echo KBPIPEOUTDIR="."      >>~/.Renviron
//     echo KBPIPEWORKDIR="./" >>~/.Renviron
//     echo RTSCRIPTS="${params.codeRdir}/" >>~/.Renviron
//
//     bash ${params.codedir}/getGridSearchStatistics.sh ${params.RTdir}/mm10 mm10
//     perl ${params.codedir}/getBestModelForEachBG.pl --s stats.mm10.tab --run
//
//     gridData=${params.accessorydir}"/RTSeq/stats.tab"
//
//     """
//     }
// }else{
//   Channel
//     .fromPath("${params.RTdir}/stats.mm10.tab")
//     .ifEmpty { exit 1, "Mouse grid stats NOT found" }
//     .set {mouseModellingStats}
//
// }

process getModelData {

  publishDir params.outdirModel,   mode: 'copy', overwrite: true
  //publishDir params.outdirRTables,   mode: 'copy', overwrite: true, pattern: '*.Rdata'

  input:
  file (af) from af17.collect()

  output:
  file("*mm10.modelMetrics.txt") into mm10ModelMetrics
  file("*.modelMetrics.txt")     into allModelMetrics
  file("*.modelMetrics.Rdata")   into allRModelMetrics

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile

  rm -f head.txt

  for stat in ${params.dataRTdir}/grid/mm10/*stats.txt ; do
    sbase=`basename \$stat`
    echo "START: \$sbase" >&2

    type=`basename \$stat |perl -pi -e 's/^(.+?)_v.+/\$1/' 2>/dev/null |perl -pi -e 's/forModel\\.mm10//' 2>/dev/null  |perl -pi -e 's/_mm10\\.forModel//' 2>/dev/null` ||true
    echo "STEP1: \$sbase" >&2

    grep -v RMSE \$stat >>\$type.mm10.modelMetrics.tmp ||true
    echo "STEP2: \$sbase" >&2

    if [ ! -f head.txt ]; then
      head -n1 \$stat >head.txt
    fi

    echo "END: \$sbase" >&2
    echo "------------------------------------------------------------------------------------" >&2
  done

  for tStat in *modelMetrics.tmp; do
    tOut=\${tStat/.tmp/.txt}
    sort -k3rn,3rn -k7rn,7rn \$tStat -o \$tStat
    cat head.txt \$tStat |perl -pi -e 's/\\s+\\/\\S+\\//\\t/g' >\$tOut
    rm \$tStat
    Rscript accessoryFiles/scripts/R/modelMetricsToRfile.R -m \$tOut
  done
  """
  }

process processBestModels {
  cache 'deep'

  tag {modelMetric}

  publishDir params.outdirModelFigs, mode: 'copy', overwrite: true, pattern: '*.png'
  publishDir params.outdirModelFigs, mode: 'copy', overwrite: true, pattern: '*WG.pdf'
  publishDir params.outdirModelFigs, mode: 'copy', overwrite: true, pattern: '*chr.pdf'
  publishDir params.outdirModelBG,   mode: 'copy', overwrite: true, pattern: '*.bedgraph'
  publishDir params.outdirModel,     mode: 'copy', overwrite: true, pattern: '*.Rdata'

  input:
  file (af) from af18.collect()
  each file(modelMetric) from allModelMetrics
  file(mm10RDBG)         from rdOrgBG_forModel_mm10.collect()
  file(oriMM)            from mouseOriginsBG_a
  file(oriMMu)           from mouseOriginsBGUnion
  // file(oriClustersPDF) from oriClustersPDF // This is just to make it wait

  output:
  file("*chr*png")                into bestModelCSPNGs
  file("*chr*pdf")                into bestModelCSPDFs
  file("*WG*png")                 into bestModelWGPNGs
  file("*WG*pdf")                 into bestModelWGPDFs
  file("*mm10.simRT.WG.bedgraph") into (bestModelSimBGs,bestModelSimBGs_a,bestModelSimBGs_b)
  file("*mm10.expRT.WG.bedgraph") into (bestModelExpBGs,bestModelExpBGs_a,bestModelExpBGs_b)
  file("*.Rdata")                 into (modelFiles, modelFiles_a, modelFiles_b, modelFiles_c)

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow ./.Rprofile
  export RTSCRIPTS="accessoryFiles/scripts/R"
  export GENOMES="accessoryFiles/genomeFiles/"

  ln -s ${params.dataRTdir}/bedgraph/final/forModel/*aph .
  ln -s ${params.annotationdir}/*bedgraph .
  ln -s accessoryFiles/ori/*omez*aph .

  bedtools shuffle -g ${mm10IDX} -i hiconf_origins.mm10.bedgraph -chrom -seed 42  |sort -k1,1 -k2n,2n -k3n,3n >randomized_hiconf_origins.mm10.bedgraph

  genome=`echo ${modelMetric} |perl -pi -e 's/^.+(hg38|mm10).+\$/\$1/'`
  name=`perl accessoryFiles/scripts/getBestModel.pl --s ${modelMetric} --g \$genome --run --bgpath . --oripath . --justgetname`

  perl accessoryFiles/scripts/getBestModel.pl --s ${modelMetric} --g \$genome --run --bgpath . --oripath . --n 100 --t 5000     >&2

  Rscript accessoryFiles/scripts/R/drawFitsFromRTmodel.R -m \$name"_"\$genome".Rdata" -d `pwd`  >&2

  """
  }

// process makeModelGIFs {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 1
//   memory '32g'
//
//   tag {model}
//
//   time { 5.hour }
//
//   module 'R/3.5.2'
//
//   publishDir params.outdirWGGIF,    mode: 'copy', overwrite: true, pattern: '*chr*.gif'
//   publishDir params.outdirCSGIF,    mode: 'copy', overwrite: true, pattern: '*WG*gif'
//
//   input:
//   each file(model) from modelFiles_a
//
//   output:
//   file('*chr*gif') into csGIFs
//   file("*WG*gif")  into wgGIFs
//
//   script:
//   """
//   ## Add globals to .Renviron (so we can see them inside R!)
//   echo KBPIPEOUTDIR="."      >>~/.Renviron
//   echo KBPIPEWORKDIR="./" >>~/.Renviron
//   echo RTSCRIPTS="${params.codeRdir}/" >>~/.Renviron
//
//   cp ${params.codeRdir}/simRT_makeGIFs.R .
//
//   Rscript simRT_makeGIFs.R -m ${model} 2>/dev/null
//
//   """
//   }
//

process getColeHiCData {

  publishDir params.outdirHiC, mode: 'copy', overwrite: true

  input:
  file (af) from af19.collect()

  output:
  file 'Zygo_HiC_100k.eigenvector.bedgraph'  into hiCZyg

  script:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122622/suppl/GSE122622_zygotene_overall.hic

  for cs in {1..19}; do
    java -Djava.awt.headless=true  -Xmx16000m  -jar \$JUICER/scripts/juicer_tools.jar eigenvector KR GSE122622_zygotene_overall.hic \$cs BP 100000 -p |perl -lane 'chomp \$_; \$_ =~ s/NaN/0/; unless (\$frm){\$frm = 1};  \$to = \$frm + 99999;  print join("\\t","chr'\$cs'",\$frm,\$to,\$_); \$frm += 100000' >>hiCZyg_100k.eigen.bedgraph
  done

  sort -k1,1 -k2n,2n -k3n,3n hiCZyg_100k.eigen.bedgraph >Zygo_HiC_100k.eigenvector.bedgraph
  """
  }

process getCASTB6hs {

  publishDir params.outdirAnnot, mode: 'copy', overwrite: true

  input:
  file (af) from af20.collect()

  output:
  file 'B6xCST.heat.bedgraph'                  into b6xcst_HeatBG
  file 'B6xCST.bias.bedgraph'                  into b6xcst_BiasBG
  file 'B6CST_B6Asymmetrics.bedgraph'         into asyHS_B6
  file 'B6CST_CSTAsymmetrics.bedgraph'        into asyHS_CST
  file 'B6CST_B6NovelAsymmetrics.bedgraph'    into asyHSNovel_B6
  file 'B6CST_CSTNovelAsymmetrics.bedgraph'   into asyHSNovel_CST

  script:
  """
  ## get B6 hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664275/suppl/GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz
  gunzip -c GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_maleHS.bedgraph

  bedtools slop -l -0.5 -r -0.5 -pct -i B6_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6_maleHS.1bp.bedgraph
  cut -f1-3 B6_maleHS.1bp.bedgraph                                                                                           >B6_maleHS.1bp.bed

  bedtools slop -l 250  -r 250       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.500bp.bedgraph
  bedtools slop -l 1500 -r 1500      -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.3Kb.bedgraph

  bedtools slop -l 500 -r 500      -i B6_maleHS.1bp.bedgraph            -g ${mm10IDX} |cut -f1-3 >B6_maleHS.1Kb.bed

  ## get CST hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1954nnn/GSM1954846/suppl/GSM1954846%5FCAST%5Fhotspots%2Etab%2Egz
  gunzip -c GSM1954846_CAST_hotspots.tab.gz |cut -f1-3,4 |grep -P \'^chr[0-9]+\' >CST_maleHS.bedgraph

  bedtools slop -l -0.5 -r -0.5 -pct -i CST_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >CST_maleHS.1bp.bedgraph
  cut -f1-3 CST_maleHS.1bp.bedgraph                                                                                        >CST_maleHS.1bp.bed

  bedtools slop -l 250  -r 250       -i CST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >CST_maleHS.500bp.bedgraph
  bedtools slop -l 1500 -r 1500      -i CST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >CST_maleHS.3Kb.bedgraph

  bedtools slop -l 500 -r 500      -i CST_maleHS.1bp.bedgraph            -g ${mm10IDX} |cut -f1-3 >CST_maleHS.1Kb.bed

  ## get CSTXb6 hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2049nnn/GSM2049312/suppl/GSM2049312%5Fdmc1hotspots%5FB6CASTF1%2EPRDM9bc%2Etxt%2Egz

  gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,\$F[2])' >B6xCST.heat.bedgraph
  gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,(\$F[3] == NA?"0.5":\$F[3]))' >B6xCST.bias.bedgraph
  gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,\$F[2],\$F[3])' >B6xCST.bg

  bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST.bg -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6xCST_maleHS.1bp.bedgraph

  bedtools slop -l 250  -r 250       -i B6xCST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6xCST_maleHS.500bp.bedgraph
  bedtools slop -l 1500 -r 1500      -i B6xCST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6xCST_maleHS.3Kb.bedgraph

  intersectBed -a B6xCST.bg -b CST_maleHS.1Kb.bed -c >C1.bg
  intersectBed -a C1.bg     -b B6_maleHS.1Kb.bed -c   >C2.bg

  cat C2.bg |perl -lane '\$type = ""; \$type = "CST" if ((\$F[5] > 0 && \$F[6] == 0) || (\$F[5] == 0 && \$F[6] == 0 && (\$F[4] > 0.75))); \$type = "B6" if ((\$F[5] == 0 && \$F[6] > 0) || (\$F[5] == 0 && \$F[6] == 0 && (\$F[4] < 0.25))); \$type = "Ambiguous" if (not \$type); print join("\\t",@F,\$type)' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_all.tab

  cat C2.bg |perl -lane '\$type = ""; \$type = "asyCST" if (\$F[4] > 0.75); \$type = "asyB6" if (\$F[4] < 0.25); if (\$type){print join("\\t",@F,\$type)}' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_allAsymmetrics.tab
  cat C2.bg |perl -lane '\$type = ""; \$type = "asyB6"  if (\$F[4] < 0.25 && \$F[5] == 0); if (\$type){print join("\\t",@F[0..3])}' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_B6Asymmetrics.bedgraph
  cat C2.bg |perl -lane '\$type = ""; \$type = "asyCST" if (\$F[4] > 0.75 && \$F[6] == 0); if (\$type){print join("\\t",@F[0..3])}' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_CSTAsymmetrics.bedgraph
  cat C2.bg |perl -lane '\$type = ""; \$type = "asyB6"  if (\$F[4] < 0.25 && \$F[5] == 0 && \$F[6] == 0); if (\$type){print join("\\t",@F[0..3])}' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_B6NovelAsymmetrics.bedgraph
  cat C2.bg |perl -lane '\$type = ""; \$type = "asyCST" if (\$F[4] > 0.75 && \$F[5] == 0 && \$F[6] == 0); if (\$type){print join("\\t",@F[0..3])}' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_CSTNovelAsymmetrics.bedgraph
  """
  }

process getHOP2hs {

  publishDir params.outdirAnnot, mode: 'copy', overwrite: true

  input:
  file (af) from af21.collect()

  output:
  file 'B6_Hop2_peaks.bedgraph' into hop2HSBG

  script:
  """
  ## get Hop2 B6 ssDNA
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3136nnn/GSM3136743/suppl/GSM3136743_Testis_SSDS_Hop2ko.ssDNA_type1.bed.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664289/suppl/GSM2664289_Testis_input_SSDS.ssDNA_type1.bed.gz

  gunzip *gz

  perl -lane '@Q = split(/_/,\$F[3]);  print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' GSM3136743_Testis_SSDS_Hop2ko.ssDNA_type1.bed >T.q30.bed
  perl -lane '@Q = split(/_/,\$F[3]);  print join("\\t",@F) if (\$Q[0] >= 30 && \$Q[1] >= 30)' GSM2664289_Testis_input_SSDS.ssDNA_type1.bed  >C.q30.bed

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 T.q30.bed |uniq >GSM3136743_Testis_SSDS_Hop2ko.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 C.q30.bed |uniq >GSM2664289_Testis_input_SSDS.bed

  n=`cat GSM3136743_Testis_SSDS_Hop2ko.bed |wc -l`
  shuf GSM2664289_Testis_input_SSDS.bed |head -n \$n |sort -k1,1 -k2n,2n -k3n,3n >GSM2664289_Testis_input_SSDS.DS.bed

  macs2 callpeak -g mm \
    -t GSM3136743_Testis_SSDS_Hop2ko.bed \\
    -c GSM2664289_Testis_input_SSDS.DS.bed \\
    --bw 1000 \\
    --keep-dup all \\
    --slocal 5000 \\
    --name B6_Hop2

  cut -f1-3 B6_Hop2_peaks.narrowPeak |grep -v ^M |grep -v chrM |sort -k1,1 -k2n,2n >B6_Hop2_peaks.bed
  perl ${params.codedir}/normalizeStrengthByAdjacentRegions.pl --bed B6_Hop2_peaks.bed --in T.q30.bed --rc --out B6_Hop2_peaks.bedgraph --tmp ./tmp
  """
  }

process makeOKSeqBW {

  publishDir params.outdirBW,   mode: 'copy', overwrite: true, pattern: '*bigwig'

  input:
  file (af)       from af38.collect()
  file(winz)      from winFiles_d.collect()

  output:
  file('*FR.bigwig')  into okSeqBW

  script:
  """
  ln -s accessoryFiles/OKseq/OKseq_ES_Petryk_SRR7535256.bam* .

  samtools view -f 16 -hb OKseq_ES_Petryk_SRR7535256.bam  >REV.bam
  samtools view -F 16 -hb OKseq_ES_Petryk_SRR7535256.bam  >FWD.bam

  bedtools bamtobed -i FWD.bam | cut -f1-3 | accessoryFiles/scripts/sortBEDByFAI.pl - ${params.genomedir}/mm10_genome.fa.fai >FWD.bed
  bedtools bamtobed -i REV.bam | cut -f1-3 | accessoryFiles/scripts/sortBEDByFAI.pl - ${params.genomedir}/mm10_genome.fa.fai >REV.bed

  mapBed -a win1ks100.bed -b FWD.bed -c 1 -o count >OkSeq.POS.frags.bg
  mapBed -a win1ks100.bed -b REV.bed -c 1 -o count >OkSeq.NEG.frags.bg

  paste OkSeq.POS.frags.bg OkSeq.NEG.frags.bg |perl -lane 'use Math::Round; \
                       \$fr = (\$F[3]+1)/(\$F[7]+1); \
                       \$logfr = log(\$fr)/log(2); \
                       \$pos=round((\$F[1]+\$F[2])/2); \
                       print join("\\t",\$F[0],\$pos-49,\$pos+50,\$logfr)' |sort -k1,1 -k2n,2n >OkSeq.FR.bg

  bedGraphToBigWig OkSeq.FR.bg ${params.genomedir}/mm10_genome.fa.fai OkSeq.FR.bigwig
  """
  }

process compareToSNSandOKSeq {

  publishDir params.outdirFigs, mode: 'copy', overwrite: true, pattern: '*png'
  publishDir params.outdirFigs, mode: 'copy', overwrite: true, pattern: '*svg'
  publishDir params.outdirFigs, mode: 'copy', overwrite: true, pattern: '*matrix*'

  input:
  file (af)       from af39.collect()
  file(oriBW)     from oriSSDSFRbw.collect()
  file mouseOriginsBG_f
  file ori_AlmeidaESC
  file ori_AlmeidaMEF
  file ori_CayrouESC
  file okSeqBW

  output:
  file('*png')     into compOriPNG
  file('*svg')     into compOriPDF
  file('*matrix*') into compOriMatrix

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow   ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"    >>./.Renviron

  sort -k1,1 -k2n,2n ${ori_CayrouESC} ${ori_AlmeidaESC} ${ori_AlmeidaMEF} >allSNS.bed

  oriESaYes=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_AlmeidaESC} -u |wc -l`
  oriESaNo=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_AlmeidaESC} -v |wc -l`
  oriMEFYes=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_AlmeidaMEF} -u |wc -l`
  oriMEFNo=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_AlmeidaMEF} -v |wc -l`
  oriEScYes=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_CayrouESC}  -u |wc -l`
  oriEScNo=`intersectBed -a ${mouseOriginsBG_f} -b ${ori_CayrouESC}  -v |wc -l`
  oriSNSYes=`intersectBed -a ${mouseOriginsBG_f} -b allSNS.bed        -u |wc -l`
  oriSNSNo=`intersectBed -a ${mouseOriginsBG_f} -b allSNS.bed        -v |wc -l`

  echo -e "sample\\toverlap\\tno"                     >oriOverlaps.txt
  echo -e "Almeida (ESC)\\t\$oriESaYes\\t\$oriESaNo" >>oriOverlaps.txt
  echo -e "Almeida (MEF)\\t\$oriMEFYes\\t\$oriMEFNo" >>oriOverlaps.txt
  echo -e "Cayrou (ESC)\\t\$oriEScYes\\t\$oriEScNo"  >>oriOverlaps.txt
  echo -e "Any\\t\$oriSNSYes\\t\$oriSNSNo"           >>oriOverlaps.txt

  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawSSDS_v_SNS_overlaps.R

  #####
  intersectBed -a ${ori_AlmeidaESC} -b ${mouseOriginsBG_f} -u |cut -f1-3 >almESCyes.bed
  intersectBed -a ${ori_AlmeidaESC} -b ${mouseOriginsBG_f} -v |cut -f1-3 >almESCno.bed

  intersectBed -a ${ori_AlmeidaMEF} -b ${mouseOriginsBG_f} -u |cut -f1-3 >almMEFyes.bed
  intersectBed -a ${ori_AlmeidaMEF} -b ${mouseOriginsBG_f} -v |cut -f1-3 >almMEFno.bed

  intersectBed -a ${ori_CayrouESC} -b ${mouseOriginsBG_f} -u |cut -f1-3 >cayESCyes.bed
  intersectBed -a ${ori_CayrouESC} -b ${mouseOriginsBG_f} -v |cut -f1-3 >cayESCno.bed

  intersectBed -a ${ori_CayrouESC} -b ${ori_AlmeidaESC} -u   |cut -f1-3 >x2ESC.bed
  intersectBed -a x2ESC.bed        -b ${mouseOriginsBG_f} -u |cut -f1-3 >x2ESCyes.bed
  intersectBed -a x2ESC.bed        -b ${mouseOriginsBG_f} -v |cut -f1-3 >x2ESCno.bed

  computeMatrix reference-point -R ${mouseOriginsBG_f} \
                                -S ${okSeqBW} mm10_WT_Rep2.FR.frags.bigwig \
                                -a 1500000 -b 1500000 \
                                --referencePoint center -p ${task.cpus} \
                                -bl accessoryFiles/blacklist/mm10.blacklist.bed \
                                -o oriSSDS_V_okSeq.matrix.gz \
                                -bs 50000 --missingDataAsZero

  plotHeatmap -m oriSSDS_V_okSeq.matrix.gz -o origins_VS_OkSeq.all.png --colorMap RdBu_r --averageTypeSummaryPlot median \
              --regionsLabel "Origins (Ori-SSDS)" \
              --samplesLabel "Ok-Seq" "Ori-SSDS" \
              --refPointLabel "0" --xAxisLabel "Distance to origin (Mb)" --yAxisLabel "Strand asymmetry; log2(F/R)" -T "" \
              --zMax 0.1 --zMin -0.1 --heatmapHeight 6

  plotHeatmap -m oriSSDS_V_okSeq.matrix.gz -o origins_VS_OkSeq.all.svg --colorMap RdBu_r --averageTypeSummaryPlot median \
              --regionsLabel "Origins (Ori-SSDS)" \
              --samplesLabel "Ok-Seq" "Ori-SSDS" \
              --refPointLabel "0" --xAxisLabel "Distance to origin (Mb)" --yAxisLabel "Strand asymmetry; log2(F/R)" -T "" \
              --zMax 0.1 --zMin -0.1 --heatmapHeight 6

  computeMatrix reference-point -R ${mouseOriginsBG_f} \
                                   almESCyes.bed \
                                   almMEFyes.bed \
                                   cayESCyes.bed\
                                   x2ESCyes.bed\
                                   almESCno.bed \
                                   almMEFno.bed \
                                   cayESCno.bed\
                                   x2ESCno.bed\
                                -S ${okSeqBW} \
                                -a 1500000 -b 1500000 \
                                --referencePoint center -p ${task.cpus} \
                                -bl accessoryFiles/blacklist/mm10.blacklist.bed \
                                -o oriVokSeq.matrix.gz \
                                -bs 50000 --missingDataAsZero

  plotHeatmap -m oriVokSeq.matrix.gz -o others_origins_VS_OkSeq.all.png --colorMap RdBu_r --averageTypeSummaryPlot median \
              --regionsLabel "Ori-SSDS" "aESC(y)" "aMEF(y)" "cESC(y)" "x2ESC(y)" "aESC(n)" "aMEF(n)" "cESC(n)" "x2ESC(n)" \
              --heatmapWidth 8 --heatmapHeight 24 --refPointLabel "0" --xAxisLabel "Distance to peak (Mb)" --yAxisLabel "Ok-Seq strand asymmetry; log2(F/R)" -T "" \
              --zMax 0.1 --zMin -0.1 --heatmapHeight 8

  plotHeatmap -m oriVokSeq.matrix.gz -o others_origins_VS_OkSeq.all.svg --colorMap RdBu_r --averageTypeSummaryPlot median \
              --regionsLabel "Ori-SSDS" "aESC(y)" "aMEF(y)" "cESC(y)" "x2ESC(y)" "aESC(n)" "aMEF(n)" "cESC(n)" "x2ESC(n)" \
              --heatmapWidth 8 --heatmapHeight 24 --refPointLabel "0" --xAxisLabel "Distance to peak (Mb)" --yAxisLabel "Ok-Seq strand asymmetry; log2(F/R)" -T "" \
              --zMax 0.1 --zMin -0.1 --heatmapHeight 8
  """
  }

process makeFigure2and4 {

  publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*.tab'
  publishDir params.outdirAnnot,   mode: 'copy', overwrite: true, pattern: '*.bedgraph'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.png'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.pdf'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.svg'

  input:
  file (af) from af22.collect()
  file(simRT)    from bestModelSimBGs.collect()
  file(expRT)    from bestModelExpBGs.collect()
  file(ori)      from mouseOriginsBG_c
  file(asyB6)    from asyHS_B6
  file(asyNov)   from asyHSNovel_B6
  file(hic)      from hiCZyg
  file(bxcSSDS)  from b6xcst_HeatBG
  file(hop2SSDS) from hop2HSBG
  file(f1PDF)    from figure1PDF

  output:
  file('*recombMetrics.bedgraph') into recombBGs
  file("rep_v_rec_MM10.tab")      into reprecTableMM10
  file("*igure2*png")             into figure2PNG
  file("*igure2*pdf")             into figure2PDF
  file("*suppFigure_allRT*png")   into figureRTSuppPNG
  file("*suppFigure_allRT*pdf")   into figureRTSuppPDF
  file("*RT_EML*png")             into figureRTEMLSuppPNG
  file("*RT_EML*pdf")             into figureRTEMLSuppPDF
  file("*igure4*png")             into figure4PNG
  file("*igure4*pdf")             into figure4PDF

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow   ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"    >>./.Renviron

  mv ${ori}   origins.recombMetrics.bedgraph
  cp accessoryFiles/recombinationData/Prdm9_AffinitySeq_AdultB6_rep1_1_peaks.bedgraph affySeq.recombMetrics.bedgraph
  cp accessoryFiles/recombinationData/Prdm9_ChIPSeq_12dppB6_rep2_1_peaks.bedgraph Prdm9ChIPSeq.recombMetrics.bedgraph
  cp accessoryFiles/recombinationData/RPA.peaks.RC.bedgraph RPA.recombMetrics.bedgraph
  mv ${asyB6}  asyB6HS.recombMetrics.bedgraph
  mv ${asyNov} asyB6Novel.recombMetrics.bedgraph
  mv ${hic} hiCZyg.recombMetrics.bedgraph
  mv ${bxcSSDS} SSDSbXc.recombMetrics.bedgraph
  mv ${hop2SSDS} SSDShop2.recombMetrics.bedgraph

  wget https://github.com/Yue-Jiang/sciliantifig/raw/master/inst/extdata/haploid.no0.hb.filtered.cast.tab

  perl -lane 'use Math::Round; \$m=round((\$F[2]+\$F[3])/2); print join("\\t",\$F[1],\$m-5000,\$m+5000,1)' haploid.no0.hb.filtered.cast.tab |\
        grep -v CHROM| \
        sort -k1,1 -k2n,2n |\
        mergeBed -i - -c 4 -o sum |\
        slopBed -i - -l 0.5 -r 0.5 -pct -g ${mm10IDX} |\
        slopBed -i - -l 1 -r 0 -g ${mm10IDX} >COYin.recombMetrics.bedgraph

  wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0492-5/MediaObjects/41586_2018_492_MOESM3_ESM.zip
  unzip 41586_2018_492_MOESM3_ESM.zip
  mv 2017-12-16654C-s3.txt BrickEtAlTable_2018.tab
  grep -vP '^(\\#|cs)' BrickEtAlTable_2018.tab |perl -lane '\$F[1] += 1499; \$F[2] -= 1500; print join("\\t",@F)' >BrickEtAlTable_2018.1bp.bed

  perl -lane 'print \$_ if (\$F[18])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,6 >SSDST1.recombMetrics.bedgraph
  perl -lane 'print \$_ if (\$F[18])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,6 >SSDST2.recombMetrics.bedgraph
  perl -lane 'print \$_ if (\$F[42])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,43 >sexBias.recombMetrics.bedgraph
  perl -lane 'print \$_ if (\$F[49])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,50 >h3k4m312DPPR1.recombMetrics.bedgraph
  perl -lane 'print \$_ if (\$F[50])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,51 >h3k4m312DPPR2.recombMetrics.bedgraph
  perl -lane 'print \$_ if (\$F[58])' BrickEtAlTable_2018.1bp.bed |cut -f1-3,55 >spo11.recombMetrics.bedgraph

  bedtools makewindows -g /data/RDCO/genomes/mm10/genome.fa.fai -w 150000 -s 50000 |sort -k1,1 -k2n,2n >ws150k50k.bed

  ## Prdm9ko
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664291/suppl/GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz
  gunzip GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz
  grep -vP 'chr[XYM]' GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph |sort -k1,1 -k2n,2n -k3n,3n >SSDSPrKO.recombMetrics.bedgraph

  ## 13R hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1954nnn/GSM1954833/suppl/GSM1954833_13R_hotspots.tab.gz
  zcat GSM1954833_13R_hotspots.tab.gz |perl -pi -e 's/\\s+(\\S)/\\t\$1/g' |sort -k1,1 -k2n,2n >SSDS13r.recombMetrics.bedgraph

  ## PRDM9h/h mouse SSDS
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2049nnn/GSM2049306/suppl/GSM2049306_dmc1hotspots_B6.PRDM9hh.txt.gz
  gunzip GSM2049306_dmc1hotspots_B6.PRDM9hh.txt.gz
  perl -lane 'print join("\\t","chr".\$F[0],\$F[1]-1500,\$F[1]+1500,\$F[2]) if (\$_ !~ /dmc1_heat/)' GSM2049306_dmc1hotspots_B6.PRDM9hh.txt |sort -k1,1 -k2n,2n -k3n,3n >SSDShh.tmp.bg
  intersectBed -a SSDShh.tmp.bg -b SSDSPrKO.recombMetrics.bedgraph -v >SSDShh.recombMetrics.bedgraph

  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1904nnn/GSM1904284/suppl/GSM1904284_H3K4me3.ForceCalledValues.B6.PRDM9hh.txt.gz
  gunzip GSM1904284_H3K4me3.ForceCalledValues.B6.PRDM9hh.txt.gz
  cut -f1-3,8 GSM1904284_H3K4me3.ForceCalledValues.B6.PRDM9hh.txt |grep -vP '^(chr\\s|\\#)' |sort -k1,1 -k2n,2n -k3n,3n >h3k4m3hh.tmp.bg
  intersectBed -a h3k4m3hh.tmp.bg -b SSDSPrKO.recombMetrics.bedgraph -v >h3k4m3hh.recombMetrics.bedgraph

  ## MYERS CO/NCOs
  wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11675-y/MediaObjects/41467_2019_11675_MOESM6_ESM.csv
  wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11675-y/MediaObjects/41467_2019_11675_MOESM7_ESM.csv

  cat 41467_2019_11675_MOESM6_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,17 |sort -k1,1 -k2n,2n |grep -v Defining            >myersCOAll.recombMetrics.bedgraph
  cat 41467_2019_11675_MOESM6_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,17 |sort -k1,1 -k2n,2n |grep -v Defining |grep HUM  >myersCOHum.recombMetrics.bedgraph
  cat 41467_2019_11675_MOESM6_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,17 |sort -k1,1 -k2n,2n |grep -v Defining |grep CAST >myersCOCAST.recombMetrics.bedgraph

  cat 41467_2019_11675_MOESM7_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,19 |sort -k1,1 -k2n,2n |grep -v allele              >myersNCOALL.recombMetrics.bedgraph
  cat 41467_2019_11675_MOESM7_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,19 |sort -k1,1 -k2n,2n |grep -v allele |grep HUM    >myersNCOHUM.recombMetrics.bedgraph
  cat 41467_2019_11675_MOESM7_ESM.csv |perl -pi -e 's/,/\\t/g' |cut -f2-4,19 |sort -k1,1 -k2n,2n |grep -v allele |grep CAST   >myersNCOCAST.recombMetrics.bedgraph

  for bg in `ls my*eco*aph`; do
    cat \$bg |perl -lane 'if (\$F[1] < \$F[2]){print join("\\t",\$F[0],\$F[1]-10,\$F[2]+10,1)}else{print join("\\t",\$F[0],\$F[2]-10,\$F[1]+10,1)}' >x.x
    sort -k1,1 -k2n,2n x.x >\$bg;
  done

  ## Get Yin et al B6xCAST B6CAST_Crossovers
  wget -O Yin_CO.B6xCAST.tab https://github.com/Yue-Jiang/sciliantifig/blob/master/inst/extdata/haploid.no0.hb.filtered.cast.tab?raw=true
  cat Yin_CO.B6xCAST.tab |perl -lane 'print join("\\t",@F[1..4],\$F[0],\$F[5])' |grep -v CHROM |sort -k1,1 -k2n,2n >Yin_CO.B6xCAST.full.bed
  cat Yin_CO.B6xCAST.tab |perl -lane '\$mid=int((\$F[2]+\$F[3])/2); print join("\\t",\$F[1],\$mid-1,\$mid,\$F[6],\$F[0],\$F[5])' |grep -v CHROM |sort -k1,1 -k2n,2n >Yin_CO.B6xCAST.midpoint.bed
  cat Yin_CO.B6xCAST.tab |perl -lane 'print join("\\t",@F[1..3],1)' |grep -v CHROM |sort -k1,1 -k2n,2n >YinCO_full.recombMetrics.bedgraph
  cat Yin_CO.B6xCAST.tab |perl -lane '\$mid=int((\$F[2]+\$F[3])/2); print join("\\t",\$F[1],\$mid-1,\$mid,1)' |grep -v CHROM |sort -k1,1 -k2n,2n >YinCO_mid.recombMetrics.bedgraph

  ## Baker B6xCAST H3K4me3 peaks
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz

  wget -O Baker_H3K4me3_BxC_peaks.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60906/suppl/GSE60906%5FH3K4me3%5FBxC%5Fmerge%5Fe%2D5%5Fpeaks%2Etxt%2Egz
  gunzip Baker_H3K4me3_BxC_peaks.txt.gz
  grep -P '^chr[0123456789]+\\s' Baker_H3K4me3_BxC_peaks.txt |cut -f1-3,6 >Baker_H3K4me3_BxC_peaks.mm9.bedgraph
  liftOver Baker_H3K4me3_BxC_peaks.mm9.bedgraph  mm9ToMm10.over.chain.gz  Baker_H3K4me3_BxC_peaks.mm10.bedgraph na
  intersectBed -a Baker_H3K4me3_BxC_peaks.mm10.bedgraph -b SSDSbXc.recombMetrics.bedgraph -wa -u >H3K4me3_BxC.recombMetrics.bedgraph

  wget -O Baker_H3K4me3_CxB_peaks.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60906/suppl/GSE60906%5FH3K4me3%5FCxB%5Fmerge%5Fe%2D5%5Fpeaks%2Etxt%2Egz
  gunzip Baker_H3K4me3_CxB_peaks.txt.gz
  grep -P '^chr[0123456789]+\\s' Baker_H3K4me3_CxB_peaks.txt |cut -f1-3,6 >Baker_H3K4me3_CxB_peaks.mm9.bedgraph
  liftOver Baker_H3K4me3_CxB_peaks.mm9.bedgraph  mm9ToMm10.over.chain.gz  Baker_H3K4me3_CxB_peaks.mm10.bedgraph na
  intersectBed -a Baker_H3K4me3_CxB_peaks.mm10.bedgraph -b SSDSbXc.recombMetrics.bedgraph -wa -u >H3K4me3_CxB.recombMetrics.bedgraph

  ## LEPTO H3K4me3 Peaks
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734408/suppl/GSM3734408%5FLE%2ER1%2EH3K4me3%2Epeaks%2Ebed%2Egz
  gunzip -c GSM3734408_LE.R1.H3K4me3.peaks.bed.gz |cut -f1-3 >Leptotene_H3K4me3.peaks.bed
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734408/suppl/GSM3734408%5FLE%2ER1%2EH3K4me3%2EmonoCorrected%2Ews25bp%2Ebigwig
  bigWigToBedGraph GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bigwig GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bg
  mapBed -sorted -a Leptotene_H3K4me3.peaks.bed -b GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bg -c 4 -o sum >h3k4m3Lep.recombMetrics.bedgraph

  ## ZYGO H3K4me3 Peaks
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414%5FZY%2ER1%2EH3K4me3%2Epeaks%2Ebed%2Egz
  gunzip -c GSM3734414_ZY.R1.H3K4me3.peaks.bed.gz |cut -f1-3 >Zygotene_H3K4me3.peaks.bed
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414%5FZY%2ER1%2EH3K4me3%2EmonoCorrected%2Ews25bp%2Ebigwig
  bigWigToBedGraph GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bigwig GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bg
  mapBed -sorted -a Zygotene_H3K4me3.peaks.bed -b GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bg -c 4 -o sum >h3k4m3Zyg.recombMetrics.bedgraph

  for rt in *mm10.expRT.WG.bedgraph;
  do
    nm="expRT_"\${rt/.expRT.WG.bedgraph/}
    ol=\${rt/.expRT.WG.bedgraph/_expRT.OL}

    perl -lane 'print \$_ if (\$F[1] >0)' \$rt >rt.bg
    sort -k1,1 -k2n,2n -k3n,3n rt.bg >\$rt
    echo \$nm >\$ol
    cut -f 4 \$rt >>\$ol

    rm -f main.bed
    echo -e "cs\tfrom\tto" >main.tab
    cut -f1-3 \$rt >>main.tab
    cut -f1-3 \$rt >>main.bed 
  done

  for rt in *mm10.simRT.WG.bedgraph;
  do
    nm="simRT_"\${rt/.simRT.WG.bedgraph/}
    ol=\${rt/.simRT.WG.bedgraph/_simRT.OL}

    perl -lane 'print \$_ if (\$F[1] >0)' \$rt >rt.bg
    sort -k1,1 -k2n,2n -k3n,3n rt.bg >\$rt
    echo \$nm >\$ol
    cut -f 4 \$rt >>\$ol

    rm -f main.bed
    echo -e "cs\tfrom\tto" >main.tab
    cut -f1-3 \$rt >>main.tab
    cut -f1-3 \$rt >>main.bed
  done

  for bgInit in *recombMetrics.bedgraph;
  do
    grep -v NA \$bgInit >ok.bg
    mv ok.bg \$bgInit

    nm=\${bgInit/.recombMetrics.bedgraph/}
    ol=\${bgInit/.recombMetrics.bedgraph/.OL}

    bedtools slop -l -0.5 -r -0.5 -pct -i \$bgInit -g ${mm10IDX} |sort -k1,1 -k2n,2n >\$nm.bgtmp
    mapBed -sorted -a main.bed -b \$nm.bgtmp -c 4 -o sum |perl -lane '\$F[3] = (\$F[3] =~ /^[\\.0123456789]+\$/ && \$F[3] ne ".")?\$F[3]:0; print join("\\t",@F)' >\$nm.bedgraph

    echo \$nm >\$ol
    cut -f 4 \$nm.bedgraph >>\$ol
  done

  #  cp ${params.codeRdir}/flipHiCdata.R .

  for bgInit in hiC*recombMetrics.bedgraph;
  do
    grep -v NA \$bgInit >ok.bg
    mv ok.bg \$bgInit

    nm=\${bgInit/.recombMetrics.bedgraph/}
    ol=\${bgInit/.recombMetrics.bedgraph/.OL}

    cp \$bgInit \$nm.bgtmp
    mapBed -sorted -a main.bed -b \$nm.bgtmp -c 4 -o sum |perl -lane '\$F[3] = (\$F[3] =~ /^[\\-\\.0123456789]+\$/ && \$F[3] ne ".")?\$F[3]:0; print join("\\t",@F)' >\$nm.bedgraph

    echo \$nm >\$ol
    cut -f 4 \$nm.bedgraph >>\$ol

    paste main.tab \$ol MeiS_EARLY_S2_2to4C_SCP3yH2AX_NA_mm10_simRT.OL |perl -pi -e 's/hiC\\S+\\tsimRT\\S+/hiC\\tRT/' |cut -f1,4,5 >hiC_v_RT.tab
    R --no-save <accessoryFiles/scripts/R/flipHiCdata.R

    paste \$nm.bedgraph \$ol |cut -f1-3,5 >tmp.bg
    mv tmp.bg \$nm.bedgraph

    mapBed -a ws150k50k.bed -b \$nm.recombMetrics.bedgraph -c 4 -o sum |grep -P '^chr[0-9]+' |perl -pi -e 's/\\t\\./\\t0/g'>\$nm.w150ks50k.bedgraph

  done

  ## Get GC content and heterochromatin content
  bedtools slop -i main.bed -g ${mm10IDX} -l -0.5 -r -0.5 -pct |bedtools slop -i - -g ${mm10IDX} -l 4999 -r 4999 >main.ok.bed
  bedtools nuc -fi ${mm10FA} -bed main.ok.bed |cut -f5 |perl -lane 'unless (\$cnt++){print "pcGC"}else{print \$_}' >GC.OL

  perl accessoryFiles/scripts/sortBEDByFAI.pl main.ok.bed ${mm10IDX} >main.faiSort.bed
  echo h3k9m3 >h3k9m3.OL
  intersectBed -a main.faiSort.bed -b ${params.annotationdir}/H3K9me3_ChIPSeq_SRR1975998.bam -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n -k3n,3n |cut -f4 >>h3k9m3.OL
  echo h3k9m2 >h3k9m2.OL
  intersectBed -a main.faiSort.bed -b ${params.annotationdir}/H3K9me2_ChIPSeq_SRR1585300.bam -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n -k3n,3n |cut -f4 >>h3k9m2.OL

  echo -e "cs\tfrom\tto" >main.tab
  cat main.ok.bed >>main.tab

  ## Add dist to P/Q tel columns
  echo -e "pDist\tqDist" >d2tels.OL
  perl accessoryFiles/scripts/d2tel.pl -i main.ok.bed -g ${mm10IDX} |cut -f4,5 >>d2tels.OL

  paste main.tab *.OL >rep_v_rec_MM10.tab

  mapBed -a ws150k50k.bed -b spo11.recombMetrics.bedgraph -c 4 -o sum |\
      grep -P '^chr[0-9]+' |perl -lane '\$F[3] = (\$F[3] eq "\\.")?0:\$F[3]; \$mid=int((\$F[2]+\$F[1])/2); print join("\\t",\$F[0],\$mid,\$F[3],"spo11",w150ks50k,"NA")' >>x.dz.tab

  cat MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10.simRT.WG.bedgraph |\
      perl -lane '\$F[3] = (\$F[3] eq "\\.")?0:\$F[3]; \$mid=int((\$F[2]+\$F[1])/2); print join("\\t",\$F[0],\$mid,\$F[3],"RT",w150ks50k,"NA")' >>x.dz.tab

  echo -e 'cs\\tpos\\tcover\\tname\\tws\\tfr' >chr12.dz.tab
  grep -wP '(cs|chr12)' x.dz.tab >>chr12.dz.tab

  #cp accessoryFiles/scripts/R/drawFigure2.R .
  #cp accessoryFiles/scripts/R/drawFNew.R .
  #cp accessoryFiles/scripts/R/drawRT_EarlyMidLate_SuppFig.R .
  cp accessoryFiles/sorting/*.fcs .
  cp accessoryFiles/img/DSBschemaBIG.png .

  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawFigure2.R >o.o 2>e.e
  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawFigure4New.R >o.o 2>e.e
  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawRT_EarlyMidLate_SuppFig.R >o.o 2>e.e
  """
  }

process makeFigure5 {

  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.png'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.pdf'

  input:
  file (af) from af27.collect()
  file(ori)      from mouseOriginsBG_f
  file(models)   from modelFiles_c.collect()

  output:
  file("*igure5*png")             into figure5PNG
  file("*igure5*pdf")             into figure5PDF

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow   ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"    >>./.Renviron

  cp accessoryFiles/blacklist/mm10.blacklist.bed .

  sort -k1,1 ${mm10IDX} |perl -lane 'print join("\\t",\$F[0],1,\$F[1])' >cs.mm10.bed

  wget https://github.com/Yue-Jiang/sciliantifig/raw/master/inst/extdata/haploid.no0.hb.filtered.cast.tab
  perl -lane 'use Math::Round; \$m=round((\$F[2]+\$F[3])/2); print join("\\t",\$F[1],\$m-5000,\$m+5000)' haploid.no0.hb.filtered.cast.tab | grep -v CHROM | sort -k1,1 -k2n,2n > B6CAST_Crossoverss.bed

  wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1954nnn/GSM1954839/suppl/GSM1954839_B6fXCASTm_hotspots.tab.gz
  gunzip -c GSM1954839_B6fXCASTm_hotspots.tab.gz | intersectBed -a - -b mm10.blacklist.bed -v > B6fXCASTm_hotspots.bedgraph

  wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz
  gunzip -c gap.txt.gz |  cut -f2-4 | sort -k1,1 -k2n,2n | mergeBed -i - | bedtools subtract -a cs.mm10.bed -b - > Gapped.mm10.bed

  bedtools nuc -fi ${mm10FA} -bed cs.mm10.bed | sed '1d'  | cut -f1,5 | sort -k1,1 -k2n,2n > cs.gc.bed

  cp MeiS_ALL_S4_2to4C_STRA8_DMRT1_mm10.Rdata rtSimModel.Rdata

  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawFigure5.R
  """
  }

process plotModelGridSearchResults {

  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.png'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*.pdf'
  publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*.csv'
  //publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*.Rdata'

  input:
  file (af) from af23.collect()
  file(rMetrics) from allRModelMetrics.collect()
  file(models)   from modelFiles_b.collect()
  //file(csGIFs)      from csGIFs

  output:
  file("Pratto*png")   into modellingMetricsPNG
  file("Pratto*pdf")   into modellingMetricsPDF
  file("summary*csv")  into ridgePlotCSVs
  //file("*Rdata") into modellingMetricsRdata

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow   ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"     >>./.Renviron

  cp accessoryFiles/sorting/*.fcs .

  R --no-save <accessoryFiles/scripts/R/drawFigure3MM10New.R

  """
  }

// Time to get some human stuff going
if (params.getHumanDMC1){
  Channel.from( ["dmc1SSDS_human_AA1" , "SRR1528821"           ,"NA"],
                ["inputSSDS_human_AA1", "SRR1528822"           ,"NA"],
                ["dmc1SSDS_human_AA2" , "SRR1528830,SRR1528831","NA"],
                ["inputSSDS_human_AA2", "SRR1528832"           ,"NA"],
                ["dmc1SSDS_human_CL4" , "NA"                   ,"dmc1SSDS_human_CL4"],
                ["inputSSDS_human_CL4", "NA"                   ,"inputSSDS_human_CL4"])
  .into {humanSSDSSRC}

  process getHumanDMC1SSDS{

      tag {name}

      publishDir params.outdir, mode: 'copy', overwrite: true

      input:
      file (af) from af24.collect()
      set(val(name),val(sra),val(fqName)) from humanSSDSSRC

      output:
      file '*e1.bed'     into humanDMC1SSDSBeds
      file '*e1.bam'     into humanDMC1SSDSBams
      file '*e1.bam.bai' into humanDMC1SSDSBais

      script:
      """
      if [ "${sra}" == "NA" ]; then

        ln -s accessoryFiles/data/hsSSDS/hg38/${fqName}.SSDS.hg38.R1.fastq .
        ln -s accessoryFiles/data/hsSSDS/hg38/${fqName}.SSDS.hg38.R2.fastq .

        nextflow run \$NXF_PIPEDIR/SSDSPipeline_1.5.groovy  \
          -c \$NXF_PIPEDIR/nextflow.local.config  \
          --fq1 ${fqName}.SSDS.hg38.R1.fastq \
          --fq2 ${fqName}.SSDS.hg38.R2.fastq \
          --r1Len 36 --r2Len 40 \
          --bwaSplitSz 20000000   \
          --name ${name}  \
          --project ${name} \
          --threads  16 \
          --sample_name ${name} \
          --library ${name} \
          --rundate  20191105 \
          --outdir . \
          --mem 32G \
          --genome hg38 \
          -with-dag nextflowDAG.dot \
          -with-report Report.html \
          -with-trace Trace.html \
          -with-timeline Timeline.html
      else
        nextflow run \$NXF_PIPEDIR/SSDSPipeline_1.5.groovy  \
          -c \$NXF_PIPEDIR/nextflow.local.config  \
          --sra ${sra} \
          --r1Len 36 --r2Len 40 \
          --bwaSplitSz 20000000   \
          --name ${name}  \
          --project ${name} \
          --threads  16 \
          --sample_name ${name} \
          --library ${name} \
          --rundate  20191105 \
          --outdir . \
          --mem 32G \
          --genome hg38 \
          -with-dag nextflowDAG.dot \
          -with-report Report.html \
          -with-trace Trace.html \
          -with-timeline Timeline.html
      fi
      """
    }
  }

Channel.from( ["AA1","AA2","CL4"] )
  .set {humanSSDSPRDM9genotypes}

if (params.callHumanHS){
  process callHotspotsHumanDMC1SSDS{

    tag {name}

    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    file (af) from af25.collect()
    val (prdm9GT) from humanSSDSPRDM9genotypes
    file (beds)   from humanDMC1SSDSBeds.collect()

    output:
    file '*peaks.bed'      into humanDMC1hotspotsBED
    file '*peaks.bedgraph' into humanDMC1hotspotsBEDGRAPH

    script:
    """
    nextflow run $NXF_PIPEDIR/callDMC1peaks.groovy \
      -c $NXF_PIPEDIR/nextflow.local.config \
      --tbed dmc1SSDS_human_${prdm9GT}.ssDNA_type1.bed \
      --cbed inputSSDS_human_${prdm9GT}.ssDNA_type1.bed \
      --genome hg38 \
      --name human_${prdm9GT} \
      --reps 3 \
      -with-report report.html \
      -with-trace trace.html \
      -with-timeline timeline.html \
      -work-dir work \
      --outdir .
    """
    }
}else{
  Channel
    .fromPath("${params.accessorydir}/recombinationData/hg38/*.peaks.RC.bedgraph")
    .ifEmpty { exit 1, "Human hotspots NOT found" }
    .set {hg38DSBhotspots}
}

Channel
  .fromPath("${params.RTdir}/bedgraph/final/forModel/*hg38*aph")
  .ifEmpty { exit 1, "RD.org bedgraphs NOT found [${params.dataRTdir}/bedgraph/final/forModel/*hg38*aph]" }
  .set {rtexp_hg38}

process makeHg38Figures {

  publishDir params.outdirRTables, mode: 'copy', overwrite: true, pattern: '*.tab'
  publishDir params.outdirAnnot,   mode: 'copy', overwrite: true, pattern: '*.bedgraph'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*igu*.png'
  publishDir params.outdirFigs,    mode: 'copy', overwrite: true, pattern: '*igu*.pdf'

  input:
  file (af) from af26.collect()
  file(expRTBG) from rtexp_hg38.collect()
  file(rdRTBG)  from rdOrgBG_forModel_hg38.collect()
  file(hs)      from hg38DSBhotspots.collect()
  file(mm10tab) from reprecTableMM10

  output:
  file('*recombMetrics.bedgraph') into humanRecombBGs
  file("rep_v_rec_HG38.tab")      into humanReprecTable

  file("*igure6*png")             into figure6PNG
  file("*igure6*pdf")             into figure6PDF

  script:
  """
  ## Get Project .Rprofile file
  cp accessoryFiles/scripts/R/Rprofile.workflow    ./.Rprofile
  echo RTSCRIPTS="accessoryFiles/scripts/R/"     >>./.Renviron

  grep -wP 'chr([\\d]*)' ${hg38IDX} >hg38Autosomes.fa.fai
  bedtools makewindows -g hg38Autosomes.fa.fai -w 10000 -s 10000 | perl -lane 'print join("\\t",\$F[0],\$F[1]+1,\$F[2])' |sort -k1,1 -k2n,2n -k3n,3n >main.bed

  echo -e "cs\tfrom\tto" >main.tab
  cat main.bed          >>main.tab

  mv human_AA1.peaks.RC.bedgraph SSDSaa1.recombMetrics.bedgraph
  mv human_AA2.peaks.RC.bedgraph SSDSaa2.recombMetrics.bedgraph
  mv human_CL4.peaks.RC.bedgraph SSDScl4.recombMetrics.bedgraph

  git clone https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination.git
  tar -zxvf ./Bherer_etal_SexualDimorphismRecombination/Refined_genetic_map_b37.tar.gz 2>e.e >o.o

  for f in ./Refined_genetic_map_b37/male_chr*txt; do
    perl -lane 'next if (\$_ =~ /pos/); \$to=\$to?\$to:0; \$from=(\$to+1); \$to = \$F[1]; print join("\\t",\$F[0],\$from,\$to,\$F[2])' \$f
  done |sort -k1,1 -k2n,2n -k3n,3n >maleRR.hg19.bedgraph

  for f in ./Refined_genetic_map_b37/female_chr*txt; do
    perl -lane 'next if (\$_ =~ /pos/); \$to=\$to?\$to:0; \$from=(\$to+1); \$to = \$F[1]; print join("\\t",\$F[0],\$from,\$to,\$F[2])' \$f
  done |sort -k1,1 -k2n,2n -k3n,3n >femaleRR.hg19.bedgraph

  wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz
  gunzip hg19ToHg38.over.chain.gz

  liftOver maleRR.hg19.bedgraph hg19ToHg38.over.chain maleRRBherer.recombMetrics.tmp unmapped
  liftOver femaleRR.hg19.bedgraph hg19ToHg38.over.chain femaleRRBherer.recombMetrics.tmp unmapped

  cat maleRRBherer.recombMetrics.tmp   |perl -lane '\$cMMb = \$F[3]; \$Mb = (\$F[2]-\$F[1])/1000000; \$cM=\$cMMb*\$Mb; print join("\\t",@F[0..2],\$cM)' >maleRRBherer.recombMetrics.bedgraph
  cat femaleRRBherer.recombMetrics.tmp |perl -lane '\$cMMb = \$F[3]; \$Mb = (\$F[2]-\$F[1])/1000000; \$cM=\$cMMb*\$Mb; print join("\\t",@F[0..2],\$cM)' >femaleRRBherer.recombMetrics.bedgraph

  wget https://science.sciencemag.org/highwire/filestream/721792/field_highwire_adjunct_files/2/aau1043_DataS1.gz
  wget https://science.sciencemag.org/highwire/filestream/721792/field_highwire_adjunct_files/3/aau1043_DataS2.gz

  zcat aau1043_DataS1.gz |grep -P '^chr\\S' |sort -k1,1 -k2n,2n |cut -f1-4 |perl -lane '\$cMMb = \$F[3]; \$Mb = (\$F[2]-\$F[1])/1000000; \$cM=\$cMMb*\$Mb; print join("\\t",@F[0..2],\$cM)' >maleRRdecode.recombMetrics.bedgraph
  zcat aau1043_DataS2.gz |grep -P '^chr\\S' |sort -k1,1 -k2n,2n |cut -f1-4 |perl -lane '\$cMMb = \$F[3]; \$Mb = (\$F[2]-\$F[1])/1000000; \$cM=\$cMMb*\$Mb; print join("\\t",@F[0..2],\$cM)' >femaleRRdecode.recombMetrics.bedgraph

  wget https://media.nature.com/original/nature-assets/ng/journal/v48/n11/extref/ng.3669-S2.xlsx

  ## Get NCO data
  ## NCOs from ChIP-based approach
  accessoryFiles/scripts/xlsx2csv.py ng.3669-S2.xlsx -s 2 -d "\\t" |cut -f3,4,5,8 |grep 1\$ |grep ^M |cut -f2-4 |perl -lane 'print join("\\t",\$F[0],\$F[1]-50,\$F[1]+50,\$F[2])' |sort -k1,1 -k2n,2n >chip_NCOs_female.bg
  accessoryFiles/scripts/xlsx2csv.py ng.3669-S2.xlsx -s 2 -d "\\t" |cut -f3,4,5,8 |grep 1\$ |grep ^P |cut -f2-4 |perl -lane 'print join("\\t",\$F[0],\$F[1]-50,\$F[1]+50,\$F[2])' |sort -k1,1 -k2n,2n >chip_NCOs_male.bg

  ## NCOs from Seq-based approach
  accessoryFiles/scripts/xlsx2csv.py ng.3669-S2.xlsx -s 3 -d "\\t" |cut -f3,4,5,8 |grep 1\$ |grep ^M |cut -f2-4 |perl -lane 'print join("\\t",\$F[0],\$F[1]-50,\$F[1]+50,\$F[2])' |sort -k1,1 -k2n,2n >seq_NCOs_female.bg
  accessoryFiles/scripts/xlsx2csv.py ng.3669-S2.xlsx -s 3 -d "\\t" |cut -f3,4,5,8 |grep 1\$ |grep ^P |cut -f2-4 |perl -lane 'print join("\\t",\$F[0],\$F[1]-50,\$F[1]+50,\$F[2])' |sort -k1,1 -k2n,2n >seq_NCOs_male.bg

  ##merge both
  sort -k1,1 -k2n,2n chip_NCOs*bg seq_NCOs*bg               >allNCOs.recombMetrics.bedgraph
  sort -k1,1 -k2n,2n chip_NCOs_female.bg seq_NCOs_female.bg >femaleNCOs.recombMetrics.bedgraph
  sort -k1,1 -k2n,2n chip_NCOs_male.bg seq_NCOs_male.bg     >maleNCOs.recombMetrics.bedgraph

  for bg in *forModel.bedgraph;
  do
    r1=\${bg/.hg38.RT.forModel.bedgraph/.RT.bedgraph}
    rt=\${r1/_hg38.forModel.bedgraph/.RT.bedgraph}
    lnk=`readlink \$bg`
    cp \$lnk \$rt

    n1="expRT_"\${rt/.RT.bedgraph/}
    nm=\${n1/expRT_RT_/expRT_}
    ol=\${rt/.RT.bedgraph/_expRT.OL}

    bedtools slop -l -0.5 -r -0.5 -pct -i \$rt -g hg38Autosomes.fa.fai |sort -k1,1 -k2n,2n -k3n,3n >\$nm.bgtmp
    mapBed -sorted -a main.bed -b \$nm.bgtmp -c 4 -o mean |perl -lane '\$F[3] = (\$F[3] =~ /^[\\-\\.0123456789]+\$/ && \$F[3] ne ".")?\$F[3]:0; print join("\\t",@F)' >\$rt

    echo \$nm >\$ol
    cut -f 4 \$rt >>\$ol
  done

  for bgInit in *recombMetrics.bedgraph;
  do
    grep -v NA \$bgInit >ok.bg
    mv ok.bg \$bgInit

    nm=\${bgInit/.recombMetrics.bedgraph/}
    ol=\${bgInit/.recombMetrics.bedgraph/.OL}

    bedtools slop -l -0.5 -r -0.5 -pct -i \$bgInit -g ${hg38IDX} |sort -k1,1 -k2n,2n >\$nm.bgtmp
    mapBed -sorted -a main.bed -b \$nm.bgtmp -c 4 -o sum |perl -lane '\$F[3] = (\$F[3] =~ /^[\\.0123456789]+\$/ && \$F[3] ne ".")?\$F[3]:0; print join("\\t",@F)' >\$nm.bedgraph

    echo \$nm >\$ol
    cut -f 4 \$nm.bedgraph >>\$ol
  done

  ## Get GC content and heterochromatin content
  bedtools slop -i main.bed -g ${hg38IDX} -l -0.5 -r -0.5 -pct |bedtools slop -i - -g ${hg38IDX} -l 4999 -r 4999 >main.ok.bed
  bedtools nuc -fi ${hg38FA} -bed main.ok.bed |cut -f5 |perl -lane 'unless (\$cnt++){print "pcGC"}else{print \$_}' >GC.OL

  perl accessoryFiles/scripts/sortBEDByFAI.pl main.ok.bed ${hg38IDX} >main.faiSort.bed

  echo -e "cs\tfrom\tto" >main.tab
  cat main.ok.bed >>main.tab

  ## Add dist to P/Q tel columns
  echo -e "pDist\tqDist" >d2tels.OL
  perl accessoryFiles/scripts/d2tel.pl -i main.ok.bed -g ${hg38IDX} |cut -f4,5 >>d2tels.OL

  paste main.tab *.OL >rep_v_rec_HG38.tab

  cat MeiS_ALL_hsS1_2to4C_SCP3_NA.RT.bedgraph |\
  perl -lane '\$F[3] = (\$F[3] eq "\\.")?0:\$F[3]; \$mid=int((\$F[2]+\$F[1])/2); print join("\\t",\$F[0],\$mid,\$F[3],"RT",w150ks50k,"NA")' >>x.dz.tab

  echo -e 'cs\\tpos\\tcover\\tname\\tws\\tfr' >chr12.dz.tab
  grep -wP '(cs|chr12)' x.dz.tab >>chr12.dz.tab

  #cp accessoryFiles/scripts/drawFigure6.R .

  R --silent --quiet --no-save <accessoryFiles/scripts/R/drawFigure6.R >&2

  """
  }
