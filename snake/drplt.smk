rule tally_distinct_observed_positions:
    output:
        "data/ucfmt.genotype.distinct_positions_tally.tsv",
    input:
        "data/ucfmt.db",
    shell:
        """
        sqlite3 -readonly -separator '\t' -header {input} > {output} <<EOF
        SELECT
          species_id,
          COUNT(DISTINCT species_position) AS tally_distinct_positions
        FROM snp_x_lib
        JOIN lib USING (lib_id)
        WHERE lib_type = 'droplet'
        GROUP BY species_id
        ;
        EOF
        """


rule tally_distinct_reference_positions:
    output:
        "data/ucfmt.reference.distinct_positions_tally.tsv",
    input:
        "data/ucfmt.db",
    shell:
        """
        sqlite3 -readonly -separator '\t' -header {input} > {output} <<EOF
        SELECT species_id, COUNT(species_position) AS tally_position
        FROM snp
        GROUP BY species_id
        ;
        EOF
        """


rule extract_barcode_to_sample:
    output:
        "data/ucfmt.barcode_to_sample.tsv",
    input:
        "data/ucfmt.db",
    shell:
        """
        sqlite3 -readonly -separator '\t' -header {input} > {output} <<EOF
        SELECT
          lib_id AS barcode,
          sample_id
        FROM lib
        ;
        EOF
        """


checkpoint extract_scg_coverage_table:
    output:
        "data/ucfmt.genotype.horizontal_coverage.tsv",
    input:
        script="scripts/extract_scg_coverages.py",
        db="data/ucfmt.db",
        aggregated="data/ucfmt.genotypes.distinct_positions_tally.tsv",
        reference="data/ucfmt.reference.distinct_positions_tally.tsv",
    params:
        threshold_horizontal_coverage=0.01,
        threshold_reference_positions=200,
    shell:
        """
        {input.script} {input.db} {input.aggregated} {input.reference} {params.threshold_horizontal_coverage} {params.threshold_reference_positions} {output}
        """


def checkpoint_extract_scg_coverage_table(wildcards, threshold):
    return (
        pd.read_table(
            checkpoints.extract_scg_coverage_table.get(**wildcards).output[0],
            index_col=["lib_id", "species_id"],
            squeeze=True,
        )[lambda x: x > threshold]
        .reset_index()
        .species_id.value_counts()
        .index.astype(str).to_list()
    )
