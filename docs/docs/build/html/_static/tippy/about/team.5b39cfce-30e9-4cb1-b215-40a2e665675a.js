selector_to_html = {"a[href=\"#team\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Team<a class=\"headerlink\" href=\"#team\" title=\"Permalink to this heading\">#</a></h1><p>The TrajAtlas project was developed in the Zhanglab at Wuhan University. We greatly appreciate every code contribution we have received thus far.\nIf you are interested in participating in this project, please check out the <a class=\"reference internal\" href=\"../contributing.html\"><span class=\"doc\">Contributing Guide</span></a> or contact us at <a class=\"reference external\" href=\"mailto:GilbertHan1011%40gmail.com\">GilbertHan1011<span>@</span>gmail<span>.</span>com</a>.</p>", "a[href=\"#core-development-team\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Core Development Team<a class=\"headerlink\" href=\"#core-development-team\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"../contributing.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Contributing Guide<a class=\"headerlink\" href=\"#contributing-guide\" title=\"Permalink to this heading\">#</a></h1><p>Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.</p><p>Please feel to reach out via opening a new <a class=\"reference external\" href=\"https://github.com/GilbertHan1011/toothAtlasManuscript/issues/new/choose\">issue</a>.</p>", "a[href=\"#principal-investigators\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Principal Investigators<a class=\"headerlink\" href=\"#principal-investigators\" title=\"Permalink to this heading\">#</a></h2>"}
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
