selector_to_html = {"a[href=\"Level1_anno_epi.html#low-quality-clusters\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Low quality clusters<a class=\"headerlink\" href=\"#low-quality-clusters\" title=\"Permalink to this heading\">#</a></h2><p>C5-2, C5-4 samples are too homogeneous. And C5-5 subcluster is too small (located in the middle of the far right side of the figure)</p><p><img alt=\"png\" src=\"../_images/C5-project.png\"/></p>", "a[href=\"Level1_anno_epi.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Level 1 Annotation<a class=\"headerlink\" href=\"#level-1-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the epithelium into two major clusters and other low quality clusters. Base on ameloblast markers, these two clusters can be annotated as ameloblasts, Non-ameloblasts.</p>", "a[href=\"#details-of-each-level-annotation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Details of Each Level Annotation<a class=\"headerlink\" href=\"#details-of-each-level-annotation\" title=\"Permalink to this heading\">#</a></h2><p>In this section, we explain the details of each level annotation.\nYou can explore the marker genes of each cluster in the <a class=\"reference external\" href=\"https://zyflab.shinyapps.io/tooth/\">shiny app</a>.</p>", "a[href=\"Level1_anno_epi.html#motivation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the epithelium into two major clusters and other low quality clusters. Base on ameloblast markers, these two clusters can be annotated as ameloblasts, Non-ameloblasts.</p>", "a[href=\"Level1_anno_epi.html#ameloblast-and-non-ameloblast\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Ameloblast and Non-ameloblast<a class=\"headerlink\" href=\"#ameloblast-and-non-ameloblast\" title=\"Permalink to this heading\">#</a></h2><p>These clusters can be devided by ameloblast markers <a class=\"reference external\" href=\"https://www.ncbi.nlm.nih.gov/gene/258\"><em>Ambn</em></a>.\n<img alt=\"png\" src=\"../_images/Epi-c5.png\"/></p>", "a[href=\"#annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Annotation<a class=\"headerlink\" href=\"#annotation\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we will adapt the annotation pipeline from <a class=\"reference external\" href=\"https://github.com/lsteuernagel/scHarmonization/\">scHarmonization</a> <span id=\"id1\">[<a class=\"reference internal\" href=\"../references.html#id22\" title=\"Lukas Steuernagel, Brian Y. H. Lam, Paul Klemm, Georgina K. C. Dowsett, Corinna A. Bauder, John A. Tadross, Tamara Sotelo Hitschfeld, Almudena del Rio Martin, Weiyi Chen, Alain J. de Solis, Henning Fenselau, Peter Davidsen, Irene Cimino, Sara N. Kohnke, Debra Rimmington, Anthony P. Coll, Andreas Beyer, Giles S. H. Yeo, and Jens C. Br\u00fcning. HypoMap\u2014a unified single-cell gene expression atlas of the murine hypothalamus. Nature Metabolism, 4(10):1402\u20131419, October 2022. doi:10.1038/s42255-022-00657-y.\">Steuernagel <em>et al.</em>, 2022</a>]</span> to our data.\nYou can just run this pipeline directly. Here, I just break down the pipeline into several steps for better understanding.</p>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(`div.content ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: false,
                onShow(instance) {MathJax.typesetPromise([instance.popper]).then(() => {});},
            });
        };
    };
    console.log("tippy tips loaded!");
};
