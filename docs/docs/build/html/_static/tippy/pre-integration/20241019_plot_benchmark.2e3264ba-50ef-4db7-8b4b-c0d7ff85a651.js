selector_to_html = {"a[href=\"#Plot-benchmark-results\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Plot benchmark results<a class=\"headerlink\" href=\"#Plot-benchmark-results\" title=\"Permalink to this heading\">#</a></h1><p>In this script, we ultilized <a class=\"reference external\" href=\"https://github.com/YosefLab/scib-metrics/tree/main\">scib-metrics</a> to bechmark and</p>"}
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
