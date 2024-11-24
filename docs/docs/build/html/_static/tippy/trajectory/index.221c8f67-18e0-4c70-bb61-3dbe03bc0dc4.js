selector_to_html = {"a[href=\"#trajectory\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory<a class=\"headerlink\" href=\"#trajectory\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we focus on odontogenesis and amelogenesis -\ntwo critical developmental processes in tooth formation.\nOur aim is to explore and understand the dynamic patterns of gene expression\nduring these processes. A key question we seek to address is whether there exists a\nuniversal molecular axis that governs mineralization across different contexts.</p>"}
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
