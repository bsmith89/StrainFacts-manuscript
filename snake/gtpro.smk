rule download_gtpro_reference_core_snps:
    output: 'raw/gtpro_refs/variation_in_species/{species_id}/core_snps.vcf.gz'
    params:
        url=lambda w: f'https://fileshare.czbiohub.org/s/waXQzQ9PRZPwTdk/download?path=%2Fvariation_in_species%2F{w.species_id}&files=core_snps.vcf.gz'
    container: None
    shell: wget_recipe

