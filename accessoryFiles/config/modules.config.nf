// PrattoEtAl config file for SLURM

process {
	//MODULE DETAILS
	withName: makeWinFiles { module = 'bedtools/2.27.1' }
	withName: getOtherOrigins { module = 'bedtools/2.27.1:ucsc/393' }
	withName: getCpGIslands { module = 'bedtools/2.27.1:ucsc/393' }
	withName: getAnnotationFiles { module = 'macs/2.1.2:bedtools/2.27.1:R/3.6.0:meme/5.0.1:ucsc/393' }
	withName: convertReplicationDomainsOrgBGs_Mouse { module = 'bedtools/2.27.1' }
	withName: convertReplicationDomainsOrgBGs_Human { module = 'bedtools/2.27.1:ucsc/393' }
	//withName: getOrigins { module = 'nextflow/20.01.0:bedtools/2.27.1:picard/2.9.2:R/3.6.0' }
	withName: getBendabilityAtOrigins { module = 'bedtools/2.27.1:samtools/1.9:R/3.6.0' }
	withName: getCoverageAtOrigins { module = 'deeptools/3.0.1:bedtools/2.27.1:samtools/1.9:ucsc/393' }
	withName: plotCoverageVG4 { module = 'nextflow/20.01.0:bedtools/2.27.1:picard/2.9.2:R/3.6.0' }
	withName: makeSlices { module = 'bedtools/2.27.1:ucsc/393:R/3.6.0' }
	withName: drawRNAhydrolysisFig { module = 'ucsc/393:R/3.6.0' }
	withName: drawSNSPeakCallingProblemsFig { module = 'ucsc/393:bedtools/2.27.1:R/3.6.0' }
	withName: originsVCpGs { module = 'bedtools/2.27.1:ucsc/393:R/3.6.0' }
	withName: analyzeOriClusters { module = 'bedtools/2.27.1:ucsc/393:R/3.6.0:deeptools/3.0.1' }
	withName: makeFigure1 { module = 'R/3.6.0:ucsc/393:bedtools/2.27.1' }
	withName: makeMouseRTSeqTable { module = 'bedtools/2.27.1:ucsc/393:R/3.6.0' }
	withName: getModelData { module = 'R/3.6.0' }
	withName: processBestModels { module = 'R/3.6.0:bedtools/2.27.1' }
	withName: getColeHiCData { module = 'juicer/1.5.6:bedtools/2.27.1:R/3.6.0:meme/5.0.1:ucsc/393' }
	withName: getCASTB6hs { module = 'macs/2.1.2:bedtools/2.27.1:R/3.6.0:meme/5.0.1:ucsc/393' }
	withName: getHOP2hs { module = 'macs/2.1.2:bedtools/2.27.1:R/3.6.0:meme/5.0.1:ucsc/393' }
	withName: makeFigure2and4 { module = 'R/3.6.0:ucsc/393:bedtools/2.27.1' }
	withName: plotModelGridSearchResults { module = 'bedtools/2.27.1:ucsc/393:R/3.6.0' }
	withName: getHumanDMC1SSDS { module = 'bedtools/2.27.1:ucsc/393' }
	withName: callHotspotsHumanDMC1SSDS { module = 'bedtools/2.27.1:ucsc/393' }
  withName: getTFBSfiles { module = 'R/3.6.0:bedtools/2.27.1:samtools/1.9' }
	withName: makeROCForTF { module = 'R/3.6.0:bedtools/2.27.1:samtools/1.9' }

	withName: makeHg38Figures { module = 'R/3.6.0:ucsc/393:bedtools/2.27.1' }
	withName:makeOKSeqBW  { module = 'bedtools/2.27.1:samtools/1.9:ucsc/393' }
  withName:compareToSNSandOKSeq { module = 'R/3.6.0:deeptools/3.0.1:bedtools/2.27.1:samtools/1.9:ucsc/393' }

	//for Call Origins workflow
	withName: shufBEDs { module = "bedtools/2.27.1:samtools/1.9" }
	withName: callReplicationOrigins { module = "bedtools/2.27.1:macs/2.1.2:samtools/1.9" }
	withName: makeOriginCallingSaturationCurve { module = "bedtools/2.27.1:R/3.6.0" }
	withName: processAndKeepOrigins { module = "bedtools/2.27.1:R/3.6.0" }
	withName: mergeOrigins { module = "bedtools/2.27.1" }
	withName: getOriginOverlaps { module = "bedtools/2.27.1:R/3.6.0" }
	withName: recalcOriginStrength { module = "bedtools/2.27.1:R/3.6.0" }
	withName:finalMergeForOrigins { module = "bedtools/2.27.1:picard/2.9.2:R/3.6.0" }
}
