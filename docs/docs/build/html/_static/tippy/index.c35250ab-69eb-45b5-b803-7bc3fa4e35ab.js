selector_to_html = {"a[href=\"#contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Contents<a class=\"headerlink\" href=\"#contents\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>", "a[href=\"#welcome-to-scatlas-s-documentation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Welcome to scAtlas\u2019s documentation!<a class=\"headerlink\" href=\"#welcome-to-scatlas-s-documentation\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>"}
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
