profiles {
  modules {
    process {
      withName: shufBEDs { module = "bedtools/2.29.0:samtools/1.9" }
      withName: callReplicationOrigins { module = "bedtools/2.29.0:macs/2.1.2:samtools/1.9" }
      withName: makeOriginCallingSaturationCurve { module = "bedtools/2.29.0:R/3.6.0" }
      withName: processAndKeepOrigins { module = "bedtools/2.29.0:R/3.6.0" }
      withName: mergeOrigins { module = "bedtools/2.29.0" }
      withName: getOriginOverlaps { module = "bedtools/2.29.0:R/3.6.0" }
      withName: recalcOriginStrength { module = "bedtools/2.29.0:R/3.6.0" }
    }
  }
}

process {
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
    time = { 1.hour * task.attempt }
  }

}
