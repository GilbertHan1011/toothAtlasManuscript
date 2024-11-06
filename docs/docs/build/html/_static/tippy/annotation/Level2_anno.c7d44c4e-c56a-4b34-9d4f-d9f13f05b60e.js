selector_to_html = {"a[href=\"#id1\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"headerlink\" href=\"#id1\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#level-2-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Level 2 annotation<a class=\"headerlink\" href=\"#level-2-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 2 annotation, we devided the mesenchyme into 9 clusters, including four large populations (\u201cC9-1\u201d,\u201dC9-2\u201d,\u201dC9-6\u201d,\u201dC9-7\u201d) and five small populations. In this level, we want to assign the consensus annotation from literature to these clusters.</p>", "a[href=\"#motivation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 2 annotation, we devided the mesenchyme into 9 clusters, including four large populations (\u201cC9-1\u201d,\u201dC9-2\u201d,\u201dC9-6\u201d,\u201dC9-7\u201d) and five small populations. In this level, we want to assign the consensus annotation from literature to these clusters.</p>"}
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
