selector_to_html = {"a[href=\"#level-1-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Level 1 Annotation<a class=\"headerlink\" href=\"#level-1-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the mesenchyme into two main groups. We are surprised to find that these two clusters are closely associated with ages.\n<img alt=\"png\" src=\"../_images/anno_c1_age.png\"/></p><p>Therefore, we named these two clusters as <code class=\"docutils literal notranslate\"><span class=\"pre\">Embryonic</span> <span class=\"pre\">Mesenchyme</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">Postnatal</span> <span class=\"pre\">Mesenchyme</span></code>.</p>", "a[href=\"#motivation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 1 annotation, we devided the mesenchyme into two main groups. We are surprised to find that these two clusters are closely associated with ages.\n<img alt=\"png\" src=\"../_images/anno_c1_age.png\"/></p><p>Therefore, we named these two clusters as <code class=\"docutils literal notranslate\"><span class=\"pre\">Embryonic</span> <span class=\"pre\">Mesenchyme</span></code> and <code class=\"docutils literal notranslate\"><span class=\"pre\">Postnatal</span> <span class=\"pre\">Mesenchyme</span></code>.</p>", "a[href=\"#postnatal-mesenchyme\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Postnatal Mesenchyme<a class=\"headerlink\" href=\"#postnatal-mesenchyme\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#embryonic-mesenchyme\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Embryonic Mesenchyme<a class=\"headerlink\" href=\"#embryonic-mesenchyme\" title=\"Permalink to this heading\">#</a></h2><p>The embryonic mesenchyme originates from cranial neural crest cells (CNCCs). After migrating into the oral region of the first pharyngeal arch, these postmigratory CNCCs commit to the dental mesenchymal lineage.</p>"}
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
