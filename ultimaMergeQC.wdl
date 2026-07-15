version 1.0

# ultimaMergeQC
#
# Merge a set of (already coordinate-sorted) Ultima CRAMs into a single sample-level
# CRAM, mark duplicates, and collect WGS QC metrics.
#
# Process:
#   1. Merge by partition, across lanes. Partitions are chr1..chr22, chrX, chrY, chrM
#      each on their own, plus one OTHER partition holding all remaining contigs and the
#      unmapped reads.
#   2. Mark duplicates within each partition (RAM scaled by the partition's size).
#   3. Merge the duplicate-marked partitions into the final sample CRAM.
#
# Ultima reads are single-end, so duplicates are defined by a read's own 5' position and
# strand (no mate). Partitioning by mapped position is therefore exact - there are no
# read pairs to split across partitions - and no rname==rnext / SPLIT routing is needed.
#
# Adapter contamination is the PCT_ADAPTER field, in the alignment_summary_metrics file 
# of collectAggregationMetrics outputs

struct CramAndIndex {
  File cram
  File crai
}

struct GenomeResources {
  String refFasta
  String refFai
  String wgsIntervalList
  String genomeModule
}

workflow ultimaMergeQC {
  input {
    Array[CramAndIndex] inputCrams
    String outputFileNamePrefix
    String reference = "hg38"
    String? ugAdapter
    String intervalsToParallelizeByString = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,OTHER"
    Float maxDuplication = 0.30
    Float maxChimerism = 0.15
    String? outputDirectory
  }

  parameter_meta {
    inputCrams: "Coordinate-sorted Ultima CRAMs (with .crai), one per lane, to merge into a single sample."
    outputFileNamePrefix: "Prefix (base file name) for all output files."
    reference: "Reference build key into the resources map (only 'hg38' is configured)."
    ugAdapter: "Optional Ultima adapter sequence; when set, passed as ADAPTER_SEQUENCE to CollectAlignmentSummaryMetrics for the PCT_ADAPTER metric. Ultima CRAMs are usually already adapter-trimmed, so residual adapter is ~0 and this can be omitted."
    intervalsToParallelizeByString: "Comma-separated partitions to scatter by. OTHER collects all non-standard contigs plus unmapped reads."
    maxDuplication: "Duplication rate above which the sample is flagged as an outlier."
    maxChimerism: "Chimerism rate (PCT_CHIMERAS) above which the sample is flagged as an outlier."
    outputDirectory: "Absolute path (on a filesystem visible from the compute nodes, e.g. /scratch/.../output) to copy the final workflow outputs into. Used as a substitute for Cromwell's final_workflow_outputs_dir."
  }

  Map[String, GenomeResources] resources = {
    "hg38": {
      "refFasta": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta",
      "refFai": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta.fai",
      "wgsIntervalList": "$HG38_ULTIMA_ROOT/wgs_calling_regions.hg38.interval_list",
      "genomeModule": "hg38-ultima/v0"
    }
  }
  GenomeResources rr = resources[reference]

  # Split the interval string into a scatterable array and compute per-partition RAM
  # coefficients (partition_size / genome_total) from the reference .fai.
  call prepareIntervals {
    input:
      str = intervalsToParallelizeByString,
      refFai = rr.refFai,
      referenceModule = rr.genomeModule
  }

  Array[String] intervals = flatten(prepareIntervals.intervals)

  scatter (interval in intervals) {
    # Step 1: subset every lane to this partition (parallel, one job per lane).
    scatter (c in inputCrams) {
      call subsetByInterval {
        input:
          inputCram = c.cram,
          inputCramIndex = c.crai,
          interval = interval,
          ncBed = prepareIntervals.ncBed,
          refFasta = rr.refFasta,
          referenceModule = rr.genomeModule,
          outputFileNamePrefix = "~{basename(c.cram, ".cram")}.~{interval}"
      }
    }

    # Step 2: merge the per-lane subsets into one partition CRAM. Skipped for a single
    # lane (nothing to merge) - that lane's subset is used directly.
    if (length(inputCrams) > 1) {
      call mergeCrams as mergeLanes {
        input:
          inputCrams = subsetByInterval.subsetCram,
          refFasta = rr.refFasta,
          referenceModule = rr.genomeModule,
          outputFileNamePrefix = "~{outputFileNamePrefix}.~{interval}.merged"
      }
    }
    File partitionCram = select_first([mergeLanes.mergedCram, subsetByInterval.subsetCram[0]])

    # Step 3: mark duplicates within the partition. RAM scales with the partition coefficient.
    call markDuplicates {
      input:
        inputCram = partitionCram,
        refFasta = rr.refFasta,
        referenceModule = rr.genomeModule,
        outputFileNamePrefix = "~{outputFileNamePrefix}.~{interval}.dupmarked",
        scaleCoefficient = prepareIntervals.intervalCoefficients[interval]
    }
  }

  # Step 4: merge the duplicate-marked partitions into one sample CRAM.
  call mergeCrams {
    input:
      inputCrams = markDuplicates.dedupCram,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      outputFileNamePrefix = outputFileNamePrefix
  }

  # QC on the merged, duplicate-marked CRAM (independent; run in parallel by the engine).
  call collectDuplicateMetrics {
    input:
      inputCram = mergeCrams.mergedCram,
      inputCramIndex = mergeCrams.mergedCramIndex,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call collectWgsMetrics {
    input:
      inputCram = mergeCrams.mergedCram,
      inputCramIndex = mergeCrams.mergedCramIndex,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      wgsIntervalList = rr.wgsIntervalList,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call collectRawWgsMetrics {
    input:
      inputCram = mergeCrams.mergedCram,
      inputCramIndex = mergeCrams.mergedCramIndex,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      wgsIntervalList = rr.wgsIntervalList,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call collectAggregationMetrics {
    input:
      inputCram = mergeCrams.mergedCram,
      inputCramIndex = mergeCrams.mergedCramIndex,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      ugAdapter = ugAdapter,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call collectReadLengthDistribution {
    input:
      inputCram = mergeCrams.mergedCram,
      inputCramIndex = mergeCrams.mergedCramIndex,
      refFasta = rr.refFasta,
      referenceModule = rr.genomeModule,
      outputFileNamePrefix = outputFileNamePrefix
  }

  Array[File] finalOutputs = [
    mergeCrams.mergedCram,
    mergeCrams.mergedCramIndex,
    collectDuplicateMetrics.metrics,
    collectWgsMetrics.metrics,
    collectRawWgsMetrics.metrics,
    collectAggregationMetrics.alignmentSummaryMetrics,
    collectAggregationMetrics.gcBiasSummaryMetrics,
    collectAggregationMetrics.gcBiasDetailMetrics,
    collectAggregationMetrics.qualityDistributionMetrics,
    collectReadLengthDistribution.readLengthDistribution,
    collectReadLengthDistribution.samtoolsStats
  ]

  if (defined(outputDirectory)) {
    call copyOutputs {
      input:
        files = finalOutputs,
        outputDirectory = outputDirectory,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }


  output {
    File mergedCram = mergeCrams.mergedCram
    File mergedCramIndex = mergeCrams.mergedCramIndex
    File duplicateMetrics = collectDuplicateMetrics.metrics
    File wgsMetrics = collectWgsMetrics.metrics
    File rawWgsMetrics = collectRawWgsMetrics.metrics
    File alignmentSummaryMetrics = collectAggregationMetrics.alignmentSummaryMetrics
    File gcBiasSummaryMetrics = collectAggregationMetrics.gcBiasSummaryMetrics
    File gcBiasDetailMetrics = collectAggregationMetrics.gcBiasDetailMetrics
    File qualityDistributionMetrics = collectAggregationMetrics.qualityDistributionMetrics
    File readLengthDistribution = collectReadLengthDistribution.readLengthDistribution
    File samtoolsStats = collectReadLengthDistribution.samtoolsStats
    File? copiedOutputsManifest = copyOutputs.copyManifest
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "Merge already-sorted Ultima CRAMs into one sample-level CRAM, mark duplicates per interval, and collect WGS QC metrics."
    dependencies: [
      {name: "gatk/4.6.2.0", url: "https://gatk.broadinstitute.org"},
      {name: "samtools/1.16.1", url: "http://www.htslib.org/"},
      {name: "rstats/4.2", url: "https://www.r-project.org/"}
    ]
    output_meta: {
      mergedCram: {description: "Merged, duplicate-marked, coordinate-sorted CRAM.", vidarr_label: "mergedCram"},
      mergedCramIndex: {description: "Index (.crai) for the merged CRAM.", vidarr_label: "mergedCramIndex"},
      duplicateMetrics: {description: "Aggregate duplicate metrics (PERCENT_DUPLICATION) over the merged CRAM.", vidarr_label: "duplicateMetrics"},
      wgsMetrics: {description: "Picard CollectWgsMetrics (Q20/MQ20-filtered coverage).", vidarr_label: "wgsMetrics"},
      rawWgsMetrics: {description: "Picard CollectRawWgsMetrics (unfiltered coverage).", vidarr_label: "rawWgsMetrics"},
      alignmentSummaryMetrics: {description: "CollectAlignmentSummaryMetrics, including PCT_ADAPTER and PCT_CHIMERAS.", vidarr_label: "alignmentSummaryMetrics"},
      gcBiasSummaryMetrics: {description: "CollectGcBiasMetrics summary.", vidarr_label: "gcBiasSummaryMetrics"},
      gcBiasDetailMetrics: {description: "CollectGcBiasMetrics detail.", vidarr_label: "gcBiasDetailMetrics"},
      qualityDistributionMetrics: {description: "QualityScoreDistribution metrics.", vidarr_label: "qualityDistributionMetrics"},
      readLengthDistribution: {description: "Read length distribution (RL section of samtools stats); Ultima fragment-size proxy.", vidarr_label: "readLengthDistribution"},
      samtoolsStats: {description: "Full samtools stats output for the merged CRAM.", vidarr_label: "samtoolsStats"}
    }
  }
}

# Split the interval string into an array and compute per-interval RAM coefficients.
task prepareIntervals {
  input {
    String str
    String refFai
    String lineSeparator = ","
    String recordSeparator = "+"
    String referenceModule
    Int jobMemory = 2
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  parameter_meta {
    str: "Interval string to split (e.g. chr1,chr2,chr3+chr4)."
    refFai: "Reference .fai (contig name in column 1, length in column 2). Read from the shared filesystem."
    lineSeparator: "Separator between intervals."
    recordSeparator: "Separator for combining multiple contigs into one interval group."
    referenceModule: "Module that provides the reference files (e.g. hg38-ultima/v0)."
    jobMemory: "Memory (GB) allocated to this job."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Tool environment modules to load (none; uses base awk/python3)."
  }

  command <<<
    set -euo pipefail

    # Intervals are separated by lineSeparator; multi-contig groups by recordSeparator.
    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t' > intervals

    # Standard contigs / keywords mentioned in the interval string (strip coordinates).
    sed 's/\t/\n/g' intervals | sed 's/:.*//' | sort -u > interval_contigs

    ### create a bed file from all contigs in the reference
    awk -v OFS="\t" '{ print $1, 1, $2 }' ~{refFai} > contigs.bed

    # OTHER bucket = every contig NOT named as a standard interval contig.
    # (More robust than grep "_" for assembly38, which also has HLA-* contigs.)
    grep -vw -F -f interval_contigs contigs.bed > other.contigs.bed || true

    python3 <<CODE
    import re
    contigs = {}
    total = 0
    with open("contigs.bed") as contigbed:
        for line in contigbed:
            contig, start, end = line.strip().split("\t")
            contigs[contig] = int(end) - int(start) + 1
            total += contigs[contig]

    with open("intervals") as interval_set, open("coefficients.txt", "w") as out:
        for line in interval_set:
            groups = line.strip().split(" ")
            interval_size = 0
            for interval in groups:
                if ":" in interval:
                    contig, start, end = re.split(r"[:-]", interval)
                    interval_size += int(end) - int(start) + 1
                elif interval == "OTHER":
                    with open("other.contigs.bed") as obed:
                        for oline in obed:
                            _, start, end = oline.strip().split("\t")
                            interval_size += int(end) - int(start) + 1
                else:
                    interval_size += contigs.get(interval, 0)
            coeff = interval_size / total if total else 0
            out.write(line.strip() + "\t" + str(coeff) + "\n")
    CODE
  >>>

  output {
    Array[Array[String]] intervals = read_tsv("intervals")
    Map[String, String] intervalCoefficients = read_map("coefficients.txt")
    File ncBed = "other.contigs.bed"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }
}

# Step 1: subset one lane to one partition. The per-lane subsets are merged implicitly
# by MarkDuplicates (multiple --INPUT), avoiding an intermediate merged CRAM.
task subsetByInterval {
  input {
    File inputCram
    File inputCramIndex
    String interval
    File ncBed
    String refFasta
    String referenceModule
    String outputFileNamePrefix
    Int jobMemory = 8
    Int cores = 4
    Int timeout = 24
    String modules = "samtools/1.16.1"
  }

  parameter_meta {
    inputCram: "One lane's input CRAM."
    inputCramIndex: "Index for inputCram."
    interval: "Chromosome name, or OTHER for all non-standard contigs plus unmapped reads."
    ncBed: "BED of non-standard contigs (used when interval == OTHER)."
    refFasta: "Reference FASTA path."
    referenceModule: "Module that provides the reference files."
    outputFileNamePrefix: "Prefix for the subset CRAM."
    jobMemory: "Memory (GB) allocated to this job."
    cores: "Threads for samtools (-@)."
    timeout: "Hours before task timeout."
    modules: "Tool environment modules to load (samtools)."
  }

  command <<<
    set -euo pipefail

    # Localize the CRAM next to its index so samtools can find the .crai.
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    if [ "~{interval}" = "OTHER" ]; then
      # All non-standard contigs plus the unmapped (no-coordinate) reads.
      samtools view -@ ~{cores} -C -T ~{refFasta} -L ~{ncBed} -o mapped.cram input.cram
      samtools view -@ ~{cores} -C -T ~{refFasta} -o unmapped.cram input.cram '*'
      samtools merge -@ ~{cores} --reference ~{refFasta} -O cram -o "~{outputFileNamePrefix}.cram" mapped.cram unmapped.cram
    else
      samtools view -@ ~{cores} -C -T ~{refFasta} -o "~{outputFileNamePrefix}.cram" input.cram "~{interval}"
    fi
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File subsetCram = "~{outputFileNamePrefix}.cram"
  }
}

# Step 2: merge the per-lane subsets (implicitly, via multiple --INPUT) and mark
# duplicates within one partition. RAM scales with the partition coefficient.
task markDuplicates {
  input {
    File inputCram
    String refFasta
    String referenceModule
    String outputFileNamePrefix
    Boolean removeDuplicates = false
    Boolean flowMode = true
    Boolean flowQIsKnownEnd = true
    Boolean flowUseUnpairedClippedEnd = true
    Boolean flowUseEndInUnpairedReads = true
    Int flowUnpairedStartUncertainty = 1
    Int flowUnpairedEndUncertainty = 0
    String? markDuplicatesAdditionalParams
    Float scaleCoefficient = 1.0
    Int jobMemory = 300
    Int minMemory = 8
    Int overhead = 4
    Int cores = 1
    Int timeout = 48
    String modules = "java/17 picard/3.4.0-patched-ultima"
  }

  parameter_meta {
    inputCram: "Merged partition CRAM (all lanes for this partition) to duplicate-mark."
    refFasta: "Reference FASTA path."
    referenceModule: "Module that provides the reference files."
    outputFileNamePrefix: "Prefix for the duplicate-marked CRAM and metrics."
    removeDuplicates: "If true, drop duplicates instead of flagging them."
    flowMode: "Ultima flow-based duplicate marking (single-end flow reads). Should stay true for Ultima data."
    flowQIsKnownEnd: "FLOW_Q_IS_KNOWN_END: treat a soft-clipped read end terminating in a quality of 0 as a known end."
    flowUseUnpairedClippedEnd: "FLOW_USE_UNPAIRED_CLIPPED_END: use the clipped, rather than aligned, end when locating unpaired reads."
    flowUseEndInUnpairedReads: "FLOW_USE_END_IN_UNPAIRED_READS: use the read end (not just start) position when keying unpaired reads."
    flowUnpairedStartUncertainty: "FLOW_UNPAIRED_START_UNCERTAINTY: positional slack (bp) allowed at the start of unpaired reads."
    flowUnpairedEndUncertainty: "FLOW_UNPAIRED_END_UNCERTAINTY: positional slack (bp) allowed at the end of unpaired reads."
    markDuplicatesAdditionalParams: "Extra arguments passed to GATK MarkDuplicates."
    scaleCoefficient: "Partition RAM scaling coefficient (partition_size / genome_total)."
    jobMemory: "Genome-wide RAM budget (GB); per-partition RAM = round(jobMemory * scaleCoefficient), floored at minMemory."
    minMemory: "Minimum RAM (GB) for any partition."
    overhead: "GB reserved for non-heap JVM overhead."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Tool environment modules to load (picard)."
  }

  Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

  command <<<
    set -euo pipefail
    # Ultima-recommended flow-based (single-end) duplicate marking.
    java -Xmx~{allocatedMemory - overhead}G -jar $PICARD_ROOT/bin/picard.jar MarkDuplicates \
      --INPUT "~{inputCram}" \
      --OUTPUT "~{outputFileNamePrefix}.cram" \
      --METRICS_FILE "~{outputFileNamePrefix}.metrics" \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --FLOW_MODE ~{flowMode} \
      --FLOW_Q_IS_KNOWN_END ~{flowQIsKnownEnd} \
      --FLOW_USE_UNPAIRED_CLIPPED_END ~{flowUseUnpairedClippedEnd} \
      --FLOW_USE_END_IN_UNPAIRED_READS ~{flowUseEndInUnpairedReads} \
      --FLOW_UNPAIRED_START_UNCERTAINTY ~{flowUnpairedStartUncertainty} \
      --FLOW_UNPAIRED_END_UNCERTAINTY ~{flowUnpairedEndUncertainty} \
      --REMOVE_DUPLICATES ~{removeDuplicates} \
      --CREATE_INDEX false \
      --VALIDATION_STRINGENCY SILENT \
      ~{markDuplicatesAdditionalParams}
  >>>

  runtime {
    memory: "~{allocatedMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File dedupCram = "~{outputFileNamePrefix}.cram"
    File metrics = "~{outputFileNamePrefix}.metrics"
  }
}

# Step 3: merge per-partition duplicate-marked CRAMs into the final sample CRAM.
task mergeCrams {
  input {
    Array[File] inputCrams
    String refFasta
    String referenceModule
    String outputFileNamePrefix
    Int jobMemory = 16
    Int cores = 8
    Int timeout = 24
    String modules = "samtools/1.16.1"
  }

  parameter_meta {
    inputCrams: "Per-partition duplicate-marked CRAMs to merge."
    refFasta: "Reference FASTA path."
    referenceModule: "Module that provides the reference files."
    outputFileNamePrefix: "Prefix for the merged CRAM."
    jobMemory: "Memory (GB) allocated to this job."
    cores: "Threads for samtools (-@)."
    timeout: "Hours before task timeout."
    modules: "Tool environment modules to load (samtools)."
  }

  command <<<
    set -euo pipefail
    samtools merge -@ ~{cores} --reference ~{refFasta} -O cram -o "~{outputFileNamePrefix}.cram" ~{sep=" " inputCrams}
    samtools index "~{outputFileNamePrefix}.cram"
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File mergedCram = "~{outputFileNamePrefix}.cram"
    File mergedCramIndex = "~{outputFileNamePrefix}.cram.crai"
  }
}

task collectDuplicateMetrics {
  input {
    File inputCram
    File inputCramIndex
    String refFasta
    String referenceModule
    String outputFileNamePrefix
    Int jobMemory = 16
    Int overhead = 4
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.6.2.0"
  }

  parameter_meta {
    inputCram: "Merged, duplicate-marked CRAM."
    inputCramIndex: "Index for inputCram."
    refFasta: "Reference FASTA path (shared filesystem)."
    referenceModule: "Module that provides the reference files."
    outputFileNamePrefix: "Prefix for the metrics file."
    jobMemory: "Memory (GB) for the job."
    overhead: "GB reserved for non-heap overhead."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Environment modules to load (gatk)."
  }

  command <<<
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    gatk --java-options "-Xmx~{jobMemory - overhead}G" CollectDuplicateMetrics \
      --INPUT input.cram \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --METRICS_FILE "~{outputFileNamePrefix}.duplicate_metrics"
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File metrics = "~{outputFileNamePrefix}.duplicate_metrics"
  }
}

task collectWgsMetrics {
  input {
    File inputCram
    File inputCramIndex
    String refFasta
    String referenceModule
    String wgsIntervalList
    String outputFileNamePrefix
    Int? readLength
    Int jobMemory = 16
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.6.2.0"
  }

  parameter_meta {
    inputCram: "Duplicate-marked CRAM."
    inputCramIndex: "Index for inputCram."
    refFasta: "Reference FASTA path (shared filesystem)."
    referenceModule: "Module that provides the reference files."
    wgsIntervalList: "Picard interval_list for the coverage territory."
    outputFileNamePrefix: "Prefix for the metrics file."
    readLength: "READ_LENGTH for the het-sensitivity calculation (default 250)."
    jobMemory: "Memory (GB) for the job."
    overhead: "GB reserved for non-heap overhead."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Environment modules to load (gatk)."
  }

  command <<<
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    gatk --java-options "-Xmx~{jobMemory - overhead}G" CollectWgsMetrics \
      --INPUT input.cram \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --OUTPUT "~{outputFileNamePrefix}.wgs_metrics.txt" \
      --INTERVALS ~{wgsIntervalList} \
      --VALIDATION_STRINGENCY SILENT \
      --INCLUDE_BQ_HISTOGRAM true \
      --USE_FAST_ALGORITHM false \
      --COUNT_UNPAIRED true \
      --COVERAGE_CAP 12500 \
      --READ_LENGTH ~{default=250 readLength}
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File metrics = "~{outputFileNamePrefix}.wgs_metrics.txt"
  }
}

task collectRawWgsMetrics {
  input {
    File inputCram
    File inputCramIndex
    String refFasta
    String referenceModule
    String wgsIntervalList
    String outputFileNamePrefix
    Int? readLength
    Int jobMemory = 16
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.6.2.0"
  }

  parameter_meta {
    inputCram: "Duplicate-marked CRAM."
    inputCramIndex: "Index for inputCram."
    refFasta: "Reference FASTA path (shared filesystem)."
    referenceModule: "Module that provides the reference files."
    wgsIntervalList: "Picard interval_list for the coverage territory."
    outputFileNamePrefix: "Prefix for the metrics file."
    readLength: "READ_LENGTH for the het-sensitivity calculation (default 250)."
    jobMemory: "Memory (GB) for the job."
    overhead: "GB reserved for non-heap overhead."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Environment modules to load (gatk)."
  }

  command <<<
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    gatk --java-options "-Xmx~{jobMemory - overhead}G" CollectRawWgsMetrics \
      --INPUT input.cram \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --OUTPUT "~{outputFileNamePrefix}.raw_wgs_metrics.txt" \
      --INTERVALS ~{wgsIntervalList} \
      --VALIDATION_STRINGENCY SILENT \
      --INCLUDE_BQ_HISTOGRAM true \
      --USE_FAST_ALGORITHM false \
      --COUNT_UNPAIRED true \
      --COVERAGE_CAP 12500 \
      --READ_LENGTH ~{default=250 readLength}
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File metrics = "~{outputFileNamePrefix}.raw_wgs_metrics.txt"
  }
}

task collectAggregationMetrics {
  input {
    File inputCram
    File inputCramIndex
    String refFasta
    String referenceModule
    String? ugAdapter
    String outputFileNamePrefix
    Int jobMemory = 16
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.6.2.0 rstats/4.2"
  }

  parameter_meta {
    inputCram: "Duplicate-marked CRAM."
    inputCramIndex: "Index for inputCram."
    refFasta: "Reference FASTA path (shared filesystem)."
    referenceModule: "Module that provides the reference files."
    ugAdapter: "Optional adapter sequence; when set, added as ADAPTER_SEQUENCE for CollectAlignmentSummaryMetrics (PCT_ADAPTER)."
    outputFileNamePrefix: "Prefix for the metrics files."
    jobMemory: "Memory (GB) for the job."
    overhead: "GB reserved for non-heap overhead."
    cores: "Cores to allocate."
    timeout: "Hours before task timeout."
    modules: "Environment modules to load (gatk, plus R for the GC-bias/quality PDFs)."
  }

  # Only pass ADAPTER_SEQUENCE when an adapter is supplied; otherwise Picard uses its default.
  String adapterArgument = if defined(ugAdapter) then "--EXTRA_ARGUMENT \"CollectAlignmentSummaryMetrics::ADAPTER_SEQUENCE=" + select_first([ugAdapter]) + "\"" else ""

  command <<<
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    gatk --java-options "-Xmx~{jobMemory - overhead}G" CollectMultipleMetrics \
      --INPUT input.cram \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --OUTPUT "~{outputFileNamePrefix}" \
      --ASSUME_SORTED true \
      --PROGRAM null \
      --PROGRAM CollectAlignmentSummaryMetrics \
      ~{adapterArgument} \
      --PROGRAM CollectGcBiasMetrics \
      --PROGRAM QualityScoreDistribution \
      --METRIC_ACCUMULATION_LEVEL SAMPLE \
      --METRIC_ACCUMULATION_LEVEL LIBRARY
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File alignmentSummaryMetrics = "~{outputFileNamePrefix}.alignment_summary_metrics"
    File gcBiasSummaryMetrics = "~{outputFileNamePrefix}.gc_bias.summary_metrics"
    File gcBiasDetailMetrics = "~{outputFileNamePrefix}.gc_bias.detail_metrics"
    File qualityDistributionMetrics = "~{outputFileNamePrefix}.quality_distribution_metrics"
  }
}

task collectReadLengthDistribution {
  input {
    File inputCram
    File inputCramIndex
    String refFasta
    String referenceModule
    String outputFileNamePrefix
    Int jobMemory = 8
    Int cores = 4
    Int timeout = 24
    String modules = "samtools/1.16.1"
  }

  parameter_meta {
    inputCram: "Duplicate-marked CRAM."
    inputCramIndex: "Index for inputCram."
    refFasta: "Reference FASTA path (shared filesystem)."
    referenceModule: "Module that provides the reference files."
    outputFileNamePrefix: "Prefix for the output files."
    jobMemory: "Memory (GB) for the job."
    cores: "Threads for samtools (-@)."
    timeout: "Hours before task timeout."
    modules: "Environment modules to load (samtools)."
  }

  command <<<
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    samtools stats -@ ~{cores} --reference ~{refFasta} input.cram > "~{outputFileNamePrefix}.samtools_stats.txt"
    # RL = read-length distribution; Ultima fragment-size proxy (GBS-7031).
    grep '^RL' "~{outputFileNamePrefix}.samtools_stats.txt" > "~{outputFileNamePrefix}.read_length_distribution.txt" || true
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{referenceModule}"
  }

  output {
    File samtoolsStats = "~{outputFileNamePrefix}.samtools_stats.txt"
    File readLengthDistribution = "~{outputFileNamePrefix}.read_length_distribution.txt"
  }
}

task copyOutputs {
  input {
    Array[File] files
    String? outputDirectory
    String outputFileNamePrefix
    Int jobMemory = 4
    Int timeout = 4
  }

  parameter_meta {
    files: "Final workflow output files to copy into outputDirectory."
    outputDirectory: "Absolute destination directory on a filesystem visible from the compute nodes. Created if it does not exist."
    outputFileNamePrefix: "Output prefix, used to name the copy manifest."
    jobMemory: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
    set -euo pipefail

    dest="~{outputDirectory}"
    mkdir -p "${dest}"

    manifest="~{outputFileNamePrefix}_copied_outputs.txt"
    : > "${manifest}"

    for f in ~{sep=' ' files}; do
      cp -f "${f}" "${dest}/"
      echo "${dest}/$(basename "${f}")" >> "${manifest}"
    done
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    File copyManifest = "~{outputFileNamePrefix}_copied_outputs.txt"
  }

  meta {
    output_meta: {
      copyManifest: "List of the destination paths the final outputs were copied to."
    }
  }
}