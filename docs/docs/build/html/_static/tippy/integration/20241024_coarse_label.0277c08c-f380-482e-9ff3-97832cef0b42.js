selector_to_html = {"a[href=\"#Validation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Validation<a class=\"headerlink\" href=\"#Validation\" title=\"Permalink to this heading\">#</a></h2><p>We ultilized marker genes of every cluster to validate our annotation results.</p>", "a[href=\"#First-level-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">First-level annotation<a class=\"headerlink\" href=\"#First-level-annotation\" title=\"Permalink to this heading\">#</a></h1><p>In this script, we aim to create level-1 annotations for downstream analysis. The level-1 annotations should align with the coarse annotations we performed earlier and distinguish the basic cell types, such as mesenchyme and epithelium</p>", "a[href=\"#Set-environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Set environment<a class=\"headerlink\" href=\"#Set-environment\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Leiden-Cluster\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Leiden Cluster<a class=\"headerlink\" href=\"#Leiden-Cluster\" title=\"Permalink to this heading\">#</a></h2><p>We gradually increased the resolution, to find the optimal resolution that can serves as the first level annotation.</p>", "a[href=\"#Save-data\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Save data<a class=\"headerlink\" href=\"#Save-data\" title=\"Permalink to this heading\">#</a></h2>"}
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
