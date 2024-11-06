selector_to_html = {"a[href=\"#why-this-matters\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Why this matters<a class=\"headerlink\" href=\"#why-this-matters\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#results\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Results<a class=\"headerlink\" href=\"#results\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#integration-methods-benchmarking\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Integration methods benchmarking<a class=\"headerlink\" href=\"#integration-methods-benchmarking\" title=\"Permalink to this heading\">#</a></h1><h2>Why this matters<a class=\"headerlink\" href=\"#why-this-matters\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#process\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Process<a class=\"headerlink\" href=\"#process\" title=\"Permalink to this heading\">#</a></h2>"}
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
