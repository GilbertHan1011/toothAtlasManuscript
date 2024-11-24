selector_to_html = {"a[href=\"#trajectory\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory<a class=\"headerlink\" href=\"#trajectory\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we focus on odontogenesis and amelogenesis -\ntwo critical developmental processes in tooth formation.\nOur aim is to explore and understand the dynamic patterns of gene expression\nduring these processes. A key question we seek to address is whether there exists a\nuniversal molecular axis that governs mineralization across different contexts.</p>", "a[href=\"multipotent.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Stemness and Start Point Inference<a class=\"headerlink\" href=\"#stemness-and-start-point-inference\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Determining the starting point of a trajectory presents a significant challenge, especially when working with integrated datasets. To address this, we can leverage key biological indicators such as cellular age and developmental potential. In this section, we will employ multiple computational methods including Cytotrace2, Cytotrace, and SCENT to systematically infer the developmental potential of cells and identify suitable trajectory starting points.</p>", "a[href=\"multipotent.html#cytotrace2\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Cytotrace2<a class=\"headerlink\" href=\"#cytotrace2\" title=\"Permalink to this heading\">#</a></h2><p><img alt=\"png\" src=\"../_images/cytotrace2.png\"/></p>", "a[href=\"multipotent.html#cytotrace\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Cytotrace<a class=\"headerlink\" href=\"#cytotrace\" title=\"Permalink to this heading\">#</a></h2><p>We applied CellRank API to infer the developmental potential of cells.\nSee this notbook for the detailed code : <a class=\"reference internal\" href=\"20241124_cytotrace.html\"><span class=\"doc std std-doc\">Estimate developmental potential inference with CytoTRACE</span></a></p>", "a[href=\"multipotent.html#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Determining the starting point of a trajectory presents a significant challenge, especially when working with integrated datasets. To address this, we can leverage key biological indicators such as cellular age and developmental potential. In this section, we will employ multiple computational methods including Cytotrace2, Cytotrace, and SCENT to systematically infer the developmental potential of cells and identify suitable trajectory starting points.</p>"}
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
