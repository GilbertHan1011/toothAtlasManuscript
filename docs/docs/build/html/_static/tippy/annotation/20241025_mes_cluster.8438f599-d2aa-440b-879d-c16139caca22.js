selector_to_html = {"a[href=\"#Leiden-Clustering-in-Mesenchyme\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Leiden Clustering in Mesenchyme<a class=\"headerlink\" href=\"#Leiden-Clustering-in-Mesenchyme\" title=\"Permalink to this heading\">#</a></h1><p>Since we subset mesenchyme from full anndata, so we rerun the <code class=\"docutils literal notranslate\"><span class=\"pre\">sc.pp.neighbors</span></code> for better clusterring.</p>"}
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
