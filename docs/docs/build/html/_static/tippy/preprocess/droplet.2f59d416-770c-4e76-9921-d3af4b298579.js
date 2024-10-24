selector_to_html = {"a[href=\"#droplet-discovery\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Droplet discovery<a class=\"headerlink\" href=\"#droplet-discovery\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>For droplet discovery, we recommended you to read these wonderful tutorials:\n<a class=\"reference external\" href=\"https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html\">Quality Control</a>;\n<a class=\"reference external\" href=\"https://bioconductor.org/books/3.19/OSCA.advanced/doublet-detection.html\">Droplet detection</a>;\nHere, we adapted scDblFinder for droplet discovery.</p>", "a[href=\"#utils-that-used-in-this-step\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Utils that used in this step<a class=\"headerlink\" href=\"#utils-that-used-in-this-step\" title=\"Permalink to this heading\">#</a></h2><p>The utils can be found in <a class=\"reference external\" href=\"https://github.com/GilbertHan1011/toothAtlasManuscript/blob/main/script/utils/seurat_utils.R\">github</a>.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>For droplet discovery, we recommended you to read these wonderful tutorials:\n<a class=\"reference external\" href=\"https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html\">Quality Control</a>;\n<a class=\"reference external\" href=\"https://bioconductor.org/books/3.19/OSCA.advanced/doublet-detection.html\">Droplet detection</a>;\nHere, we adapted scDblFinder for droplet discovery.</p>", "a[href=\"#process-scripts\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Process Scripts<a class=\"headerlink\" href=\"#process-scripts\" title=\"Permalink to this heading\">#</a></h2><p>The utils can be found in <a class=\"reference external\" href=\"https://github.com/GilbertHan1011/toothAtlasManuscript/blob/master/anno_base/20241008_droplet.R\">github</a>.</p>"}
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
