selector_to_html = {"a[href=\"#step-4-prune-tree\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Step 4: Prune Tree<a class=\"headerlink\" href=\"#step-4-prune-tree\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#script\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Script<a class=\"headerlink\" href=\"#script\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#input\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Input<a class=\"headerlink\" href=\"#input\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#output\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Output<a class=\"headerlink\" href=\"#output\" title=\"Permalink to this heading\">#</a></h2>"}
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
