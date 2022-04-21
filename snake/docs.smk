rule sort_bib_from_raw:
    output:
        "doc/bibliography.bib",
    input:
        script="scripts/sort_bib.py",
        bib=["doc/bibliography_raw.bib"],
    shell:
        "{input.script} {input.bib} > {output}"


rule render_figure_to_png:
    output:
        "fig/{stem}_figure.w{width}.png",
    input:
        "doc/static/{stem}_figure.svg",
    params:
        width=lambda w: int(w.width),
    shell:
        """
        inkscape {input} --export-width={params.width} --export-filename {output}
        """

ruleorder: render_figure_to_png > render_pdf_to_png_imagemagick


rule render_pdf_to_png_imagemagick:
    output:
        "fig/{stem}.dpi{dpi}.png",
    input:
        "fig/{stem}.pdf",
    params:
        dpi=lambda w: int(w.dpi),
    shell:
        """
        convert -units PixelsPerInch -density {params.dpi} {input} {output}
        """


rule render_pdf_to_tiff_imagemagick:
    output:
        "fig/{stem}.dpi{dpi}.tiff",
    input:
        "fig/{stem}.pdf",
    params:
        dpi=lambda w: int(w.dpi),
    shell:
        """
        convert -units PixelsPerInch -density {params.dpi} {input} {output}
        """


rule link_static_pdf_figure:
    output:
        "fig/{stem}_figure.pdf",
    input:
        "doc/static/{stem}_figure.pdf",
    shell:
        alias_recipe


rule pdf_to_eps:
    output:
        "fig/{stem}.eps",
    input:
        "fig/{stem}.pdf",
    shell:
        """
        cd fig
        pdf2ps {wildcards.stem}.pdf
        ps2eps {wildcards.stem}.ps
        rm {wildcards.stem}.ps
        """


rule render_figure_to_pdf:
    output:
        "fig/{stem}_figure.pdf",
    input:
        "doc/static/{stem}_figure.svg",
    shell:
        """
        inkscape {input} --export-filename {output}
        """


rule build_submission_docx:
    output:
        "build/{stem}.docx",
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        template="doc/static/example_style.docx",
        csl="doc/citestyle.csl",
        figures=lambda w: config["figures"][w.stem],
    conda: 'conda/pandoc.yaml'
    shell:
        """
        pandoc --from markdown --to docx \
               --standalone --self-contained --reference-doc {input.template} \
               --filter pandoc-crossref --csl {input.csl} --citeproc \
               --bibliography={input.bib} -s {input.source} -o {output}
        """


localrules:
    build_submission_docx,


rule build_submission_pdf:
    output:
        "build/{stem}.pdf",
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        csl="doc/citestyle.csl",
        figures=lambda w: config["figures"][w.stem],
    shell:
        """
        pandoc --from markdown --to pdf \
               --filter pandoc-crossref --csl {input.csl} --citeproc \
               --pdf-engine=xelatex \
               --bibliography={input.bib} -s {input.source} -o {output}
        """


localrules:
    build_submission_pdf,


rule compile_submission_folder:
    output:
        directory("build/submission"),
    input:
        docx='build/submission.docx',
        fig1_tiff='fig/compute_profiling_figure.dpi500.tiff',
        fig2_tiff='fig/benchmarks_figure.dpi500.tiff',
        fig3_tiff='fig/scg_comparison_figure.dpi500.tiff',
        fig4_tiff='fig/coclustering_figure.dpi500.tiff',
        fig5_tiff='fig/biogeography_figure.dpi500.tiff',
        fig6_tiff='fig/ld_decay_figure.dpi500.tiff',
        fig1_pdf='fig/compute_profiling_figure.pdf',
        fig2_pdf='fig/benchmarks_figure.pdf',
        fig3_pdf='fig/scg_comparison_figure.pdf',
        fig4_pdf='fig/coclustering_figure.pdf',
        fig5_pdf='fig/biogeography_figure.pdf',
        fig6_pdf='fig/ld_decay_figure.pdf',
        fig1_eps='fig/compute_profiling_figure.eps',
        fig2_eps='fig/benchmarks_figure.eps',
        fig3_eps='fig/scg_comparison_figure.eps',
        fig4_eps='fig/coclustering_figure.eps',
        fig5_eps='fig/biogeography_figure.eps',
        fig6_eps='fig/ld_decay_figure.eps',
    shell:
        """
        mkdir -p {output}
        cp {input.fig1_tiff} {output}/Figure1.tiff
        cp {input.fig2_tiff} {output}/Figure2.tiff
        cp {input.fig3_tiff} {output}/Figure3.tiff
        cp {input.fig4_tiff} {output}/Figure4.tiff
        cp {input.fig5_tiff} {output}/Figure5.tiff
        cp {input.fig6_tiff} {output}/Figure6.tiff
        cp {input.fig1_pdf} {output}/Figure1.pdf
        cp {input.fig2_pdf} {output}/Figure2.pdf
        cp {input.fig3_pdf} {output}/Figure3.pdf
        cp {input.fig4_pdf} {output}/Figure4.pdf
        cp {input.fig5_pdf} {output}/Figure5.pdf
        cp {input.fig6_pdf} {output}/Figure6.pdf
        cp {input.fig1_eps} {output}/Figure1.eps
        cp {input.fig2_eps} {output}/Figure2.eps
        cp {input.fig3_eps} {output}/Figure3.eps
        cp {input.fig4_eps} {output}/Figure4.eps
        cp {input.fig5_eps} {output}/Figure5.eps
        cp {input.fig6_eps} {output}/Figure6.eps
        cp {input.docx} {output}/Submission.docx
        """
