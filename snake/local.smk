rule import_bibliography:
    output: 'doc/bibliography_raw.bib'
    input: '/Users/byronsmith/Documents/Literature/strainfacts.bib'
    shell: alias_recipe
