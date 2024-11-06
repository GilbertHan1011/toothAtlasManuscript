selector_to_html = {"a[href=\"#first-round-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">First round annotation<a class=\"headerlink\" href=\"#first-round-annotation\" title=\"Permalink to this heading\">#</a></h1><p>In the first round annotation, we recommend to coarsely annotate the cell types. It serves two purposes:</p>"}
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
