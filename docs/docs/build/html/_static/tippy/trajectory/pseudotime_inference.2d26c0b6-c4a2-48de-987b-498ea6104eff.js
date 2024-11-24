selector_to_html = {"a[href=\"#pseudotime-inference\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Pseudotime Inference<a class=\"headerlink\" href=\"#pseudotime-inference\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Based on our previous analyses, we have identified two main developmental trajectories in tooth formation: one from bud mesenchyme to odontoblasts during embryonic development, and another from apical papilla to odontoblasts during postnatal and adult stages. To quantitatively characterize these developmental paths, we will employ Slingshot to infer pseudotime orderings along these trajectories.</p>", "a[href=\"#slingshot\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Slingshot<a class=\"headerlink\" href=\"#slingshot\" title=\"Permalink to this heading\">#</a></h2><p>We use diffusion map as input for slingshot, as we can adjust the dimension of diffusion map to tune the results to align with biological insights.\nAnd because the embryonic trajectory and postnatal trajectory have greate difference, we will apply slingshot to these two datasets separately.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Based on our previous analyses, we have identified two main developmental trajectories in tooth formation: one from bud mesenchyme to odontoblasts during embryonic development, and another from apical papilla to odontoblasts during postnatal and adult stages. To quantitatively characterize these developmental paths, we will employ Slingshot to infer pseudotime orderings along these trajectories.</p>"}
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
