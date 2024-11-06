selector_to_html = {"a[href=\"#first-round-preprocess\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">First round preprocess<a class=\"headerlink\" href=\"#first-round-preprocess\" title=\"Permalink to this heading\">#</a></h1><h2>Aim of first round preprocess<a class=\"headerlink\" href=\"#aim-of-first-round-preprocess\" title=\"Permalink to this heading\">#</a></h2><p>The goal of first round preprocess is to generate a clean and standard data for downstream analysis. So we build an standard pipeline to perform basic quality control and annotation.</p>", "a[href=\"#aim-of-first-round-preprocess\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Aim of first round preprocess<a class=\"headerlink\" href=\"#aim-of-first-round-preprocess\" title=\"Permalink to this heading\">#</a></h2><p>The goal of first round preprocess is to generate a clean and standard data for downstream analysis. So we build an standard pipeline to perform basic quality control and annotation.</p>", "a[href=\"#utils-that-used-in-this-step\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Utils that used in this step<a class=\"headerlink\" href=\"#utils-that-used-in-this-step\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#basic-process-script\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">basic process script<a class=\"headerlink\" href=\"#basic-process-script\" title=\"Permalink to this heading\">#</a></h2>"}
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
