selector_to_html = {"a[href=\"#save-data\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Save data<a class=\"headerlink\" href=\"#save-data\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#leiden-cluster\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Leiden Cluster<a class=\"headerlink\" href=\"#leiden-cluster\" title=\"Permalink to this heading\">#</a></h2><p>We gradually increased the resolution, to find the optimal resolution that can serves as the first level annotation.</p>", "a[href=\"#cluster-the-datasets-to-get-first-level-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Cluster the datasets to get first level annotation<a class=\"headerlink\" href=\"#cluster-the-datasets-to-get-first-level-annotation\" title=\"Permalink to this heading\">#</a></h1><p>In this script, we aim to create level-1 annotations for downstream analysis. The level-1 annotations should align with the coarse annotations we performed earlier and distinguish the basic cell types, such as mesenchyme and epithelium</p>", "a[href=\"#validation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Validation<a class=\"headerlink\" href=\"#validation\" title=\"Permalink to this heading\">#</a></h2><p>We ultilized marker genes of every cluster to validate our annotation results.</p>", "a[href=\"#set-environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Set environment<a class=\"headerlink\" href=\"#set-environment\" title=\"Permalink to this heading\">#</a></h2>"}
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
