selector_to_html = {"a[href=\"#Hypertune:-grid-hypertune\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Hypertune: grid hypertune<a class=\"headerlink\" href=\"#Hypertune:-grid-hypertune\" title=\"Permalink to this heading\">#</a></h1><p>The number of highly variable gene (HVG) and latent number is very important for scANVI algrithm. Therefore, we using <a class=\"reference external\" href=\"https://www.dremio.com/wiki/grid-search/\">grid search</a> For hyper-parameter tuning.</p>"}
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
