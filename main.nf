nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW }      from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM }     from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC }               from './modules/nf-core/trimmomatic/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN }  from './modules/nf-core/kraken2/kraken2/main.nf'
include { MULTIQC }                   from './modules/nf-core/multiqc/main.nf'
include { paramsSummaryMap          } from 'plugin/nf-schema'

def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params
        .keySet()
        .each { group ->
            def group_params = summary_params.get(group)
            // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>${group}</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                group_params
                    .keySet()
                    .sort()
                    .each { param ->
                        summary_section += "        <dt>${param}</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                    }
                summary_section += "    </dl>\n"
            }
        }

    def yaml_file_text = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n" as String
    yaml_file_text     += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text     += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text     += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text     += "plot_type: 'html'\n"
    yaml_file_text     += "data: |\n"
    yaml_file_text     += "${summary_section}"

    return yaml_file_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

process CHECK_OR_DOWNLOAD_DB {
    tag 'kraken_db'
    container 'docker.io/library/python:3.11-slim'
    publishDir "refs/kraken2", mode: 'copy'
    cache 'lenient'

    output:
    path "k2_minusb_20250402", emit: db_dir

    script:
    """
    set -euo pipefail
    DB_DIR='k2_minusb_20250402'
    if [ ! -d "\$DB_DIR" ]; then
        pip install --quiet --no-cache-dir --user gdown
        gdown --id 1C4aisqMEmUiIv-jNBzxTkKX1AqFKgYen -O kraken_db.tar.gz
        mkdir -p "\$DB_DIR" && tar -xzf kraken_db.tar.gz -C "\$DB_DIR"
        rm kraken_db.tar.gz
    fi
    """
}

workflow {
    // ────────────────────────────────────────────────────────────────────────────
    // 1) Setup
    // ────────────────────────────────────────────────────────────────────────────
    println "Project directory = $projectDir"
    println "Timestamp         = ${params.timestamp}"
    println "Output directory  = ${params.outdir}"

    def kraken2_db = file('refs/kraken2/k2_minusb_20250402')
    if( ! kraken2_db.exists() ) {
        CHECK_OR_DOWNLOAD_DB()
    }

    // ────────────────────────────────────────────────────────────────────────────
    // 2) Raw FASTQ discovery
    // ────────────────────────────────────────────────────────────────────────────
    Channel
      .fromFilePairs("raw_data/*_{1,2}.fastq.gz", size: 2)
      .ifEmpty { error("No paired FASTQ files found…") }
      .map { id, reads -> tuple([ id: id, single_end: false ], reads) }
      .set { raw_reads }

    raw_reads.view { meta, reads -> "Sample: ${meta.id}, files: ${reads}" }

    // ────────────────────────────────────────────────────────────────────────────
    // 3) Run modules
    // ────────────────────────────────────────────────────────────────────────────
    FASTQC_RAW(   raw_reads )
    TRIMMOMATIC(  raw_reads )
    FASTQC_TRIM(
      TRIMMOMATIC.out.trimmed_reads.map { meta, reads ->
          def nm = meta.clone(); nm.id = "${meta.id}_trimmed"
          tuple(nm, reads)
      }
    )
    KRAKEN(
      TRIMMOMATIC.out.trimmed_reads.map { meta, reads -> tuple(meta, reads) },
      kraken2_db, false, true
    )

    // ────────────────────────────────────────────────────────────────────────────
    // 4) nf-core MultiQC wiring
    // ────────────────────────────────────────────────────────────────────────────

    // 4a) config, logo, custom methods
    ch_multiqc_config = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    // 4b) workflow summary for MultiQC
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    // 4c) collate module outputs
    ch_collated_versions = Channel.from( file('versions.yml') )  // nf-core versions file

    // 4d) initial merge of QC files
    ch_multiqc_files = Channel
      .merge(
        FASTQC_RAW .out.html,
        FASTQC_RAW .out.zip,
        FASTQC_TRIM.out.html,
        FASTQC_TRIM.out.zip,
        TRIMMOMATIC .out.summary,
        KRAKEN      .out.report
      )
      .map { meta, file -> file }

    // 4e) mix in workflow summary, versions, methods
    ch_multiqc_files = ch_multiqc_files.mix(
      ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_methods_description = Channel.value(
      methodsDescriptionText(
        params.multiqc_methods_description ? file(params.multiqc_methods_description) :
        file("$projectDir/assets/methods_description_template.yml")
      )
    )
    ch_multiqc_files = ch_multiqc_files.mix(
      ch_collated_versions,
      ch_methods_description.collectFile(
        name: 'methods_description_mqc.yaml',
        sort: true
      )
    )

    // 4f) finally, run the MultiQC module
    MULTIQC(
      ch_multiqc_files.collect(),           // all files as LIST<Path>
      ch_multiqc_config.toList(),           // default config
      ch_multiqc_custom_config.toList(),    // user-supplied config
      ch_multiqc_logo.toList(),             // logo, if any
      [], []                                // replace_names & sample_names (none)
    )

    ch_multiqc_report_list = MULTIQC.out.report.toList()

    // ────────────────────────────────────────────────────────────────────────────
    // 5) Done
    // ────────────────────────────────────────────────────────────────────────────
    println("Workflow completed at: ${new Date()}")
    println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    println("Output directory: ${params.outdir}")
}
