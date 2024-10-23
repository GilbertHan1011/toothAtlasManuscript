selector_to_html = {"a[href=\"#citing-scatlas\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Citing scAtlas<a class=\"headerlink\" href=\"#citing-scatlas\" title=\"Permalink to this heading\">\u00b6</a></h1><p><a class=\"reference external\" href=\"https://www.biorxiv.org/content/10.1101/2024.05.28.596174v1.full\">Our work</a> has been posted on bioRxiv. If scAtlas is useful for your research, please consider citing our work as follows:</p>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(` ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: false,

            });
        };
    };
    console.log("tippy tips loaded!");
};
