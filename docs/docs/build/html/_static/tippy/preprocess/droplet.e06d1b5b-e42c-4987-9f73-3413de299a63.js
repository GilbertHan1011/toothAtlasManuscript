selector_to_html = {"a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>For droplet discovery, we recommended you to read these wonderful tutorials:\n<a class=\"reference external\" href=\"https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html\">Quality Control</a>;\n<a class=\"reference external\" href=\"https://bioconductor.org/books/3.19/OSCA.advanced/doublet-detection.html\">Droplet detection</a>;</p>", "a[href=\"#scripts\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Scripts<a class=\"headerlink\" href=\"#scripts\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#droplet-discovery\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Droplet discovery<a class=\"headerlink\" href=\"#droplet-discovery\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>For droplet discovery, we recommended you to read these wonderful tutorials:\n<a class=\"reference external\" href=\"https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html\">Quality Control</a>;\n<a class=\"reference external\" href=\"https://bioconductor.org/books/3.19/OSCA.advanced/doublet-detection.html\">Droplet detection</a>;</p>"}
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
