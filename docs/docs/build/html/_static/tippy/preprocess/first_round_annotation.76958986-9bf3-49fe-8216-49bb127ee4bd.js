selector_to_html = {"a[href=\"#curated-marker-genes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Curated marker genes<a class=\"headerlink\" href=\"#curated-marker-genes\" title=\"Permalink to this heading\">#</a></h2><p>Though the cell type may have great variance across tissues and ages, the main cell types are consistent across datasets. In the first round annotation, the annotation should be unified. Base on the literature, we summarized the cell types that commonly exist in the literature, and make a marker gene lists.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In the first round annotation, we recommend to coarsely annotate the cell types. It serves two purposes:</p>", "a[href=\"#first-round-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">First round annotation<a class=\"headerlink\" href=\"#first-round-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In the first round annotation, we recommend to coarsely annotate the cell types. It serves two purposes:</p>", "a[href=\"#procedure\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Procedure<a class=\"headerlink\" href=\"#procedure\" title=\"Permalink to this heading\">#</a></h2><p>Load Environment</p>"}
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
