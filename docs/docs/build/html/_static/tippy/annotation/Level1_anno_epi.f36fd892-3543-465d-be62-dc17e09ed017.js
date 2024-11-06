selector_to_html = {"a[href=\"#low-quality-clusters\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Low quality clusters<a class=\"headerlink\" href=\"#low-quality-clusters\" title=\"Permalink to this heading\">#</a></h2><p>C5-2, C5-4 samples are too homogeneous. And C5-5 subcluster is too small (located in the middle of the far right side of the figure)</p><p><img alt=\"png\" src=\"../_images/C5-project.png\"/></p>", "a[href=\"#motivation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the epithelium into two major clusters and other low quality clusters. Base on ameloblast markers, these two clusters can be annotated as ameloblasts, Non-ameloblasts.</p>", "a[href=\"#ameloblast-and-non-ameloblast\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Ameloblast and Non-ameloblast<a class=\"headerlink\" href=\"#ameloblast-and-non-ameloblast\" title=\"Permalink to this heading\">#</a></h2><p>These clusters can be devided by ameloblast markers <a class=\"reference external\" href=\"https://www.ncbi.nlm.nih.gov/gene/258\"><em>Ambn</em></a>.\n<img alt=\"png\" src=\"img/annotation_epi/Epi-C5.png\"/></p>", "a[href=\"#level-1-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Level 1 Annotation<a class=\"headerlink\" href=\"#level-1-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the epithelium into two major clusters and other low quality clusters. Base on ameloblast markers, these two clusters can be annotated as ameloblasts, Non-ameloblasts.</p>"}
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
