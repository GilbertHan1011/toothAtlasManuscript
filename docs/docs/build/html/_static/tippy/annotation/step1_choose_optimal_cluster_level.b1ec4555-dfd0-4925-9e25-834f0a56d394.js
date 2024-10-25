selector_to_html = {"a[href=\"#step-1-choose-the-optimal-cluster-level\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Step 1: Choose the optimal cluster level<a class=\"headerlink\" href=\"#step-1-choose-the-optimal-cluster-level\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In this step, we will select the optimal cluster level from the previous <a class=\"reference internal\" href=\"20241025_mes_cluster.html\"><span class=\"std std-doc\">clustering results</span></a> to annotate the mesenchyme cells.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>In this step, we will select the optimal cluster level from the previous <a class=\"reference internal\" href=\"20241025_mes_cluster.html\"><span class=\"std std-doc\">clustering results</span></a> to annotate the mesenchyme cells.</p>", "a[href=\"#process\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Process<a class=\"headerlink\" href=\"#process\" title=\"Permalink to this heading\">#</a></h2><p><img alt=\"png\" src=\"../_images/annotation_step1.png\"/></p>"}
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
