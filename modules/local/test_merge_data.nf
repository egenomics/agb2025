nextflow.enable.dsl = 2

include { HOMOSAPINENS_CONTAMINATION } from './homo_contam.nf'

workflow {
    def run_id = params.run_id ?: "TEST_RUN"

    // Define the input files
    def metadata_path = "runs/${run_id}/metadata/metadata.tsv"
    def kraken_reports_path = "runs/${run_id}/kraken/*.kraken2.report.txt"

    Channel.fromPath(metadata_path, checkIfExists: true)
        .set { metadata_ch }

    Channel.fromPath(kraken_reports_path, checkIfExists: true)
        .set { kraken_reports_ch }

    // Run only this module
    HOMOSAPINENS_CONTAMINATION(metadata_ch, kraken_reports_ch)
}
