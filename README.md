# ultimaMergeQC

Merge already-sorted Ultima CRAMs into one sample-level CRAM, mark duplicates per interval, and collect WGS QC metrics.

## Overview

## Dependencies

* [gatk 4.6.2.0](https://gatk.broadinstitute.org)
* [samtools 1.16.1](http://www.htslib.org/)
* [rstats 4.2](https://www.r-project.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run ultimaMergeQC.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputCrams`|Array[CramAndIndex]|Coordinate-sorted Ultima CRAMs (with .crai), one per lane, to merge into a single sample.
`outputFileNamePrefix`|String|Prefix (base file name) for all output files.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`reference`|String|"hg38"|Reference build key into the resources map (only 'hg38' is configured).
`ugAdapter`|String?|None|Optional Ultima adapter sequence; when set, passed as ADAPTER_SEQUENCE to CollectAlignmentSummaryMetrics for the PCT_ADAPTER metric. Ultima CRAMs are usually already adapter-trimmed, so residual adapter is ~0 and this can be omitted.
`intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,OTHER"|Comma-separated partitions to scatter by. OTHER collects all non-standard contigs plus unmapped reads.
`maxDuplication`|Float|0.3|Duplication rate above which the sample is flagged as an outlier.
`maxChimerism`|Float|0.15|Chimerism rate (PCT_CHIMERAS) above which the sample is flagged as an outlier.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`prepareIntervals.lineSeparator`|String|","|Separator between intervals.
`prepareIntervals.recordSeparator`|String|"+"|Separator for combining multiple contigs into one interval group.
`prepareIntervals.jobMemory`|Int|2|Memory (GB) allocated to this job.
`prepareIntervals.cores`|Int|1|Cores to allocate.
`prepareIntervals.timeout`|Int|1|Hours before task timeout.
`prepareIntervals.modules`|String|""|Tool environment modules to load (none; uses base awk/python3).
`subsetByInterval.jobMemory`|Int|8|Memory (GB) allocated to this job.
`subsetByInterval.cores`|Int|4|Threads for samtools (-@).
`subsetByInterval.timeout`|Int|24|Hours before task timeout.
`subsetByInterval.modules`|String|"samtools/1.16.1"|Tool environment modules to load (samtools).
`mergeLanes.jobMemory`|Int|16|Memory (GB) allocated to this job.
`mergeLanes.cores`|Int|8|Threads for samtools (-@).
`mergeLanes.timeout`|Int|24|Hours before task timeout.
`mergeLanes.modules`|String|"samtools/1.16.1"|Tool environment modules to load (samtools).
`markDuplicates.removeDuplicates`|Boolean|false|If true, drop duplicates instead of flagging them.
`markDuplicates.flowMode`|Boolean|true|Ultima flow-based duplicate marking (single-end flow reads). Should stay true for Ultima data.
`markDuplicates.flowQIsKnownEnd`|Boolean|true|FLOW_Q_IS_KNOWN_END: treat a soft-clipped read end terminating in a quality of 0 as a known end.
`markDuplicates.flowUseUnpairedClippedEnd`|Boolean|true|FLOW_USE_UNPAIRED_CLIPPED_END: use the clipped, rather than aligned, end when locating unpaired reads.
`markDuplicates.flowUseEndInUnpairedReads`|Boolean|true|FLOW_USE_END_IN_UNPAIRED_READS: use the read end (not just start) position when keying unpaired reads.
`markDuplicates.flowUnpairedStartUncertainty`|Int|1|FLOW_UNPAIRED_START_UNCERTAINTY: positional slack (bp) allowed at the start of unpaired reads.
`markDuplicates.flowUnpairedEndUncertainty`|Int|0|FLOW_UNPAIRED_END_UNCERTAINTY: positional slack (bp) allowed at the end of unpaired reads.
`markDuplicates.markDuplicatesAdditionalParams`|String?|None|Extra arguments passed to GATK MarkDuplicates.
`markDuplicates.jobMemory`|Int|300|Genome-wide RAM budget (GB); per-partition RAM = round(jobMemory * scaleCoefficient), floored at minMemory.
`markDuplicates.minMemory`|Int|8|Minimum RAM (GB) for any partition.
`markDuplicates.overhead`|Int|4|GB reserved for non-heap JVM overhead.
`markDuplicates.cores`|Int|1|Cores to allocate.
`markDuplicates.timeout`|Int|48|Hours before task timeout.
`markDuplicates.modules`|String|"gatk/4.6.2.0"|Tool environment modules to load (gatk).
`mergeCrams.jobMemory`|Int|16|Memory (GB) allocated to this job.
`mergeCrams.cores`|Int|8|Threads for samtools (-@).
`mergeCrams.timeout`|Int|24|Hours before task timeout.
`mergeCrams.modules`|String|"samtools/1.16.1"|Tool environment modules to load (samtools).
`collectDuplicateMetrics.jobMemory`|Int|16|Memory (GB) for the job.
`collectDuplicateMetrics.overhead`|Int|4|GB reserved for non-heap overhead.
`collectDuplicateMetrics.cores`|Int|1|Cores to allocate.
`collectDuplicateMetrics.timeout`|Int|24|Hours before task timeout.
`collectDuplicateMetrics.modules`|String|"gatk/4.6.2.0"|Environment modules to load (gatk).
`collectWgsMetrics.readLength`|Int?|None|READ_LENGTH for the het-sensitivity calculation (default 250).
`collectWgsMetrics.jobMemory`|Int|16|Memory (GB) for the job.
`collectWgsMetrics.overhead`|Int|6|GB reserved for non-heap overhead.
`collectWgsMetrics.cores`|Int|1|Cores to allocate.
`collectWgsMetrics.timeout`|Int|24|Hours before task timeout.
`collectWgsMetrics.modules`|String|"gatk/4.6.2.0"|Environment modules to load (gatk).
`collectRawWgsMetrics.readLength`|Int?|None|READ_LENGTH for the het-sensitivity calculation (default 250).
`collectRawWgsMetrics.jobMemory`|Int|16|Memory (GB) for the job.
`collectRawWgsMetrics.overhead`|Int|6|GB reserved for non-heap overhead.
`collectRawWgsMetrics.cores`|Int|1|Cores to allocate.
`collectRawWgsMetrics.timeout`|Int|24|Hours before task timeout.
`collectRawWgsMetrics.modules`|String|"gatk/4.6.2.0"|Environment modules to load (gatk).
`collectAggregationMetrics.jobMemory`|Int|16|Memory (GB) for the job.
`collectAggregationMetrics.overhead`|Int|6|GB reserved for non-heap overhead.
`collectAggregationMetrics.cores`|Int|1|Cores to allocate.
`collectAggregationMetrics.timeout`|Int|24|Hours before task timeout.
`collectAggregationMetrics.modules`|String|"gatk/4.6.2.0 rstats/4.2"|Environment modules to load (gatk, plus R for the GC-bias/quality PDFs).
`collectReadLengthDistribution.jobMemory`|Int|8|Memory (GB) for the job.
`collectReadLengthDistribution.cores`|Int|4|Threads for samtools (-@).
`collectReadLengthDistribution.timeout`|Int|24|Hours before task timeout.
`collectReadLengthDistribution.modules`|String|"samtools/1.16.1"|Environment modules to load (samtools).
`checkPreValidation.jobMemory`|Int|4|Memory (GB) for the job.
`checkPreValidation.cores`|Int|1|Cores to allocate.
`checkPreValidation.timeout`|Int|2|Hours before task timeout.
`checkPreValidation.modules`|String|""|Environment modules to load (none; uses system python3).


### Outputs

Output | Type | Description | Labels
---|---|---|---
`mergedCram`|File|Merged, duplicate-marked, coordinate-sorted CRAM.|vidarr_label: mergedCram
`mergedCramIndex`|File|Index (.crai) for the merged CRAM.|vidarr_label: mergedCramIndex
`duplicateMetrics`|File|Aggregate duplicate metrics (PERCENT_DUPLICATION) over the merged CRAM.|vidarr_label: duplicateMetrics
`perIntervalDuplicateMetrics`|Array[File]|Per-interval MarkDuplicates metrics files.|vidarr_label: perIntervalDuplicateMetrics
`wgsMetrics`|File|Picard CollectWgsMetrics (Q20/MQ20-filtered coverage).|vidarr_label: wgsMetrics
`rawWgsMetrics`|File|Picard CollectRawWgsMetrics (unfiltered coverage).|vidarr_label: rawWgsMetrics
`alignmentSummaryMetrics`|File|CollectAlignmentSummaryMetrics, including PCT_ADAPTER and PCT_CHIMERAS.|vidarr_label: alignmentSummaryMetrics
`gcBiasSummaryMetrics`|File|CollectGcBiasMetrics summary.|vidarr_label: gcBiasSummaryMetrics
`gcBiasDetailMetrics`|File|CollectGcBiasMetrics detail.|vidarr_label: gcBiasDetailMetrics
`qualityDistributionMetrics`|File|QualityScoreDistribution metrics.|vidarr_label: qualityDistributionMetrics
`readLengthDistribution`|File|Read length distribution (RL section of samtools stats); Ultima fragment-size proxy.|vidarr_label: readLengthDistribution
`samtoolsStats`|File|Full samtools stats output for the merged CRAM.|vidarr_label: samtoolsStats
`duplicationRate`|Float|PERCENT_DUPLICATION parsed from the duplicate metrics.|vidarr_label: duplicationRate
`chimerismRate`|Float|PCT_CHIMERAS parsed from the alignment summary metrics.|vidarr_label: chimerismRate
`isOutlierData`|Boolean|True if duplication or chimerism exceed the configured thresholds.|vidarr_label: isOutlierData


## Commands
This section lists command(s) run by ultimaMergeQC workflow

* Running ultimaMergeQC

```
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
```
```
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
```
```
    set -euo pipefail
    # Ultima-recommended flow-based (single-end) duplicate marking.
    gatk --java-options "-Xmx~{allocatedMemory - overhead}G" MarkDuplicates \
      --INPUT="~{inputCram}" \
      --OUTPUT="~{outputFileNamePrefix}.cram" \
      --METRICS_FILE="~{outputFileNamePrefix}.metrics" \
      --REFERENCE_SEQUENCE="~{refFasta}" \
      --FLOW_MODE ~{flowMode} \
      --FLOW_Q_IS_KNOWN_END ~{flowQIsKnownEnd} \
      --FLOW_USE_UNPAIRED_CLIPPED_END ~{flowUseUnpairedClippedEnd} \
      --FLOW_USE_END_IN_UNPAIRED_READS ~{flowUseEndInUnpairedReads} \
      --FLOW_UNPAIRED_START_UNCERTAINTY ~{flowUnpairedStartUncertainty} \
      --FLOW_UNPAIRED_END_UNCERTAINTY ~{flowUnpairedEndUncertainty} \
      --REMOVE_DUPLICATES=~{removeDuplicates} \
      --CREATE_INDEX=false \
      --VALIDATION_STRINGENCY=SILENT \
      ~{markDuplicatesAdditionalParams}
```
```
    set -euo pipefail
    samtools merge -@ ~{cores} --reference ~{refFasta} -O cram -o "~{outputFileNamePrefix}.cram" ~{sep=" " inputCrams}
    samtools index "~{outputFileNamePrefix}.cram"
```
```
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    gatk --java-options "-Xmx~{jobMemory - overhead}G" CollectDuplicateMetrics \
      --INPUT input.cram \
      --REFERENCE_SEQUENCE ~{refFasta} \
      --METRICS_FILE "~{outputFileNamePrefix}.duplicate_metrics"
```
```
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
```
```
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
```
```
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
```
```
    set -euo pipefail
    ln -s ~{inputCram} input.cram
    ln -s ~{inputCramIndex} input.cram.crai

    samtools stats -@ ~{cores} --reference ~{refFasta} input.cram > "~{outputFileNamePrefix}.samtools_stats.txt"
    # RL = read-length distribution; Ultima fragment-size proxy (GBS-7031).
    grep '^RL' "~{outputFileNamePrefix}.samtools_stats.txt" > "~{outputFileNamePrefix}.read_length_distribution.txt" || true
```
```
    set -euo pipefail

    grep -A 1 PERCENT_DUPLICATION ~{duplicationMetrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerismMetrics} | grep -v OF_PAIR > chimerism.csv

    python3 <<CODE
    import csv
    with open('duplication.csv') as dupfile:
        reader = csv.DictReader(dupfile, delimiter='\t')
        for row in reader:
            with open("duplication_value.txt", "w") as f:
                f.write(row['PERCENT_DUPLICATION'])

    with open('chimerism.csv') as chimfile:
        reader = csv.DictReader(chimfile, delimiter='\t')
        for row in reader:
            with open("chimerism_value.txt", "w") as f:
                f.write(row['PCT_CHIMERAS'])
    CODE
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
