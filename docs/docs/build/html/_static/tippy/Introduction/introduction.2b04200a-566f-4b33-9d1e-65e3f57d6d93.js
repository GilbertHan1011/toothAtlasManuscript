selector_to_html = {"a[href=\"#introduction-to-scatlas\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction to scAtlas<a class=\"headerlink\" href=\"#introduction-to-scatlas\" title=\"Permalink to this heading\">#</a></h1><h2>Why atlas<a class=\"headerlink\" href=\"#why-atlas\" title=\"Permalink to this heading\">#</a></h2><p>One of the most exciting techenical innovation is the rapid development of single-cell technech, enabling us to have a close look in single cell resolution.</p>", "a[href=\"#why-atlas\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Why atlas<a class=\"headerlink\" href=\"#why-atlas\" title=\"Permalink to this heading\">#</a></h2><p>One of the most exciting techenical innovation is the rapid development of single-cell technech, enabling us to have a close look in single cell resolution.</p>"}
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
