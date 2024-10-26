selector_to_html = {"a[href=\"step2_mtree.html#script\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Script<a class=\"headerlink\" href=\"#script\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"step2_mtree.html#input\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Input<a class=\"headerlink\" href=\"#input\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"step2_mtree.html#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In this step, we will use the MRTree algorithm  <span id=\"id1\">[<a class=\"reference internal\" href=\"../references.html#id23\" title=\"Minshi Peng, Brie Wamsley, Andrew G Elkins, Daniel H Geschwind, Yuting Wei, and Kathryn Roeder. Cell type hierarchy reconstruction via reconciliation of multi-resolution cluster tree. Nucleic Acids Research, 49(16):e91, June 2021. doi:10.1093/nar/gkab481.\">Peng <em>et al.</em>, 2021</a>]</span> to build a hierarchical clustering tree from the input labels. MRTree is a powerful tool using multiple resolution clustering results to build a hierarchical tree.</p>", "a[href=\"step2_mtree.html#output\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Output<a class=\"headerlink\" href=\"#output\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"step2_mtree.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Step 2: MRTree<a class=\"headerlink\" href=\"#step-2-mrtree\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In this step, we will use the MRTree algorithm  <span id=\"id1\">[<a class=\"reference internal\" href=\"../references.html#id23\" title=\"Minshi Peng, Brie Wamsley, Andrew G Elkins, Daniel H Geschwind, Yuting Wei, and Kathryn Roeder. Cell type hierarchy reconstruction via reconciliation of multi-resolution cluster tree. Nucleic Acids Research, 49(16):e91, June 2021. doi:10.1093/nar/gkab481.\">Peng <em>et al.</em>, 2021</a>]</span> to build a hierarchical clustering tree from the input labels. MRTree is a powerful tool using multiple resolution clustering results to build a hierarchical tree.</p>", "a[href=\"#annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Annotation<a class=\"headerlink\" href=\"#annotation\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we will adapt the annotation pipeline from <a class=\"reference external\" href=\"https://github.com/lsteuernagel/scHarmonization/\">scHarmonization</a> <span id=\"id1\">[<a class=\"reference internal\" href=\"../references.html#id22\" title=\"Lukas Steuernagel, Brian Y. H. Lam, Paul Klemm, Georgina K. C. Dowsett, Corinna A. Bauder, John A. Tadross, Tamara Sotelo Hitschfeld, Almudena del Rio Martin, Weiyi Chen, Alain J. de Solis, Henning Fenselau, Peter Davidsen, Irene Cimino, Sara N. Kohnke, Debra Rimmington, Anthony P. Coll, Andreas Beyer, Giles S. H. Yeo, and Jens C. Br\u00fcning. HypoMap\u2014a unified single-cell gene expression atlas of the murine hypothalamus. Nature Metabolism, 4(10):1402\u20131419, October 2022. doi:10.1038/s42255-022-00657-y.\">Steuernagel <em>et al.</em>, 2022</a>]</span> to our data.\nYou can just run this pipeline directly. Here, I just break down the pipeline into several steps for better understanding.</p>"}
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
