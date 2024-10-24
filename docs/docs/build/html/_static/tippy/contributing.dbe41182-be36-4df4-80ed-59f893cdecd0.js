selector_to_html = {"a[href=\"#contributing-guide\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Contributing Guide<a class=\"headerlink\" href=\"#contributing-guide\" title=\"Permalink to this heading\">#</a></h1><p>Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.</p><p>Please feel to reach out via opening a new <a class=\"reference external\" href=\"https://github.com/GilbertHan1011/toothAtlasManuscript/issues/new/choose\">issue</a>.</p>"}
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
