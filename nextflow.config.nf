// Nextflow config file for Pratto et al
executor = "slurm"

profiles {

  modules {
    includeConfig "modules.config.nf"
  }

  local {
    executor = "local"
  }

  singularity {
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY'
    process.container = "${params.projectdir}/accessoryFiles/singularityImg/PrattoEtAl.sif"
  }

  none {
    // Add custom configs here
  }

  callOrigins {
    includeConfig "callOris.config.nf"
  }
}

params.outdir = './output_PrattoEtAl'
timeline {
  enabled = true
  file = "${params.outdir}/PrattoEtAl_timeline.html"
}
trace {
  enabled = true
  file = "${params.outdir}/PrattoEtAl_trace.txt"
}

manifest {
  description = 'Pratto et. al analysis pipeline. Author: Kevin Brick.'
}

// PrattoEtAl config file for SLURM
process {
	project = 'PrattoEtAl_Pipeline'

	// Defaults maxs
	max_memory = 128.GB
	max_cpus = 16
	max_time = 240.h

	//DEFAULT PROCESS PROPS
	cpus = { 1 * task.attempt }
	memory = { 8.GB * task.attempt }
	time = { 2.h * task.attempt }

	errorStrategy = { task.exitStatus == 143 ? 'retry' : 'finish' }
	maxRetries = 2
	maxErrors = '-1'
	scratch = '/lscratch/$SLURM_JOBID'
	clusterOptions = ' --gres=lscratch:300 '

	//PROCESS-SPECIFIC RESOURCES
	withName:makeWinFiles {
		cpus = { 2 * task.attempt }
		memory = { 4.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:getOtherOrigins {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:getCpGIslands {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:getAnnotationFiles {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 2.hour * task.attempt }
	}
	withName:convertReplicationDomainsOrgBGs_Mouse {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 0.5.hour * task.attempt }
	}
	withName:convertReplicationDomainsOrgBGs_Human {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 0.5.hour * task.attempt }
	}
	withName:getOrigins {
		cpus = { 8 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 8.hour * task.attempt }
	}
	withName:getBendabilityAtOrigins {
		cpus = { 16 * task.attempt }
		memory = { 16.GB * task.attempt }
		time = { 2.hour * task.attempt }
	}
	withName:getCoverageAtOrigins {
		cpus = { 16 * task.attempt }
		memory = { 16.GB * task.attempt }
		time = { 2.hour * task.attempt }
	}
	withName:plotCoverageVG4 {
		cpus = { 16 * task.attempt }
		memory = { 16.GB * task.attempt }
		time = { 8.hour * task.attempt }
	}
	withName:makeSlices {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:drawRNAhydrolysisFig {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 0.2.hour * task.attempt }
	}
	withName:drawSNSPeakCallingProblemsFig {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 0.2.hour * task.attempt }
	}
	withName:originsVCpGs {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:analyzeOriClusters {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 3.hour * task.attempt }
	}
	withName:makeFigure1 {
		cpus = { 1 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 0.5.hour * task.attempt }
	}
	withName:makeMouseRTSeqTable {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:getModelGridSearchResults {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 2.h * task.attempt }
	}
	withName:getModelData {
		cpus = { 2 * task.attempt }
		memory = { 16.GB * task.attempt }
		time = { 4.hour * task.attempt }
	}
	withName:processBestModels {
		cpus = { 1 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 1.5.hour * task.attempt }
	}
	withName:makeModelGIFs {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 2.h * task.attempt }
	}
	withName:getColeHiCData {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 2.hour * task.attempt }
	}
	withName:getCASTB6hs {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 0.5.hour * task.attempt }
	}
	withName:getHOP2hs {
		cpus = { 1 * task.attempt }
		memory = { 8.GB * task.attempt }
		time = { 2.hour * task.attempt }
	}
	withName:makeFigure2and4 {
		cpus = { 1 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 2.hour * task.attempt }
    time = { 0.5.hour * task.attempt }
	}
	withName:plotModelGridSearchResults {
		cpus = { 1 * task.attempt }
		memory = { 128.GB * task.attempt }
		time = { 1.hour * task.attempt }
	}
	withName:getHumanDMC1SSDS {
		cpus = { 16 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 24.hour * task.attempt }
	}
	withName:callHotspotsHumanDMC1SSDS {
		cpus = { 16 * task.attempt }
		memory = { 32.GB * task.attempt }
		time = { 24.hour * task.attempt }
	}
	withName:makeHg38Figures {
		cpus = { 1 * task.attempt }
		memory = { 32.GB * task.attempt }
    time = { 1.hour * task.attempt }
	}

  withName:makeOKSeqBW {
		cpus = { 4 * task.attempt }
		memory = { 16.GB * task.attempt }
    time = { 2.hour * task.attempt }
	}

  withName:compareToSNSandOKSeq {
		cpus = { 16 * task.attempt }
		memory = { 32.GB * task.attempt }
    time = { 2.hour * task.attempt }
	}

  //For call Origins workflow
  withName:shufBEDs {
    cpus = { 2 * task.attempt }
    memory = { 16.GB * task.attempt }
    time = { 2.hour * task.attempt }
  }

  withName:callReplicationOrigins {
    cpus = { 4 * task.attempt }
    memory = { 32.GB * task.attempt }
    time = { 2.hour * task.attempt }
  }

  withName:makeOriginCallingSaturationCurve {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
  }

  withName:processAndKeepOrigins {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 2.hour * task.attempt }
  }

  withName:mergeOrigins {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
  }

  withName:getOriginOverlaps {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
  }

  withName:recalcOriginStrength {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
  }

  withName:finalMergeForOrigins {
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }
    time = { 2.hour * task.attempt }
  }
}
