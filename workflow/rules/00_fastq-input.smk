# ───────────────────────────────────────────────
# Rule: Stage FASTQ (intermediates)
# ───────────────────────────────────────────────
rule stage_flnc_part:
    message: "Stage FASTQ: {wildcards.sample}, part {wildcards.part}"
    input:
        src = lambda wc: PART_MAP[(wc.sample, int(wc.part))]
    output:
        dest = temp("01_mapping/{sample}/flnc_parts/flnc_part{part}.fastq.gz")
    benchmark:
        "benchmarks/01_mapping/{sample}_stage_fastq_{part}.txt"
    log:
        "logs/01_mapping/{sample}_stage_fastq_{part}.log"
    threads: 1
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        r'''
        (
          echo "Copying/gzipping FASTQ part to staging"
          case "{input.src}" in
            *.fastq.gz)  cp -f "{input.src}" "{output.dest}" ;;
            *.fastq)     gzip -c "{input.src}" > "{output.dest}" ;;
            *)           echo "Unsupported FASTQ extension" >&2; exit 2 ;;
          esac
          echo "Done"
        ) &> "{log}"
        '''
