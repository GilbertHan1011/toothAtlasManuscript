selector_to_html = {"a[href=\"pseudotime_inference.html#slingshot\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Slingshot<a class=\"headerlink\" href=\"#slingshot\" title=\"Permalink to this heading\">#</a></h2><p>We use diffusion map as input for slingshot, as we can adjust the dimension of diffusion map to tune the results to align with biological insights.\nAnd because the embryonic trajectory and postnatal trajectory have greate difference, we will apply slingshot to these two datasets separately.</p>", "a[href=\"#trajectory\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory<a class=\"headerlink\" href=\"#trajectory\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we focus on odontogenesis and amelogenesis -\ntwo critical developmental processes in tooth formation.\nOur aim is to explore and understand the dynamic patterns of gene expression\nduring these processes. A key question we seek to address is whether there exists a\nuniversal molecular axis that governs mineralization across different contexts.</p>", "a[href=\"pseudotime_inference.html#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Based on our previous analyses, we have identified two main developmental trajectories in tooth formation: one from bud mesenchyme to odontoblasts during embryonic development, and another from apical papilla to odontoblasts during postnatal and adult stages. To quantitatively characterize these developmental paths, we will employ Slingshot to infer pseudotime orderings along these trajectories.</p>", "a[href=\"pseudotime_inference.html#result\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Result<a class=\"headerlink\" href=\"#result\" title=\"Permalink to this heading\">#</a></h2><p><img alt=\"png\" src=\"../_images/slingshot_merge.png\"/>\n<img alt=\"png\" src=\"../_images/slingshot_split.png\"/></p>", "a[href=\"pseudotime_inference.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Pseudotime Inference<a class=\"headerlink\" href=\"#pseudotime-inference\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Based on our previous analyses, we have identified two main developmental trajectories in tooth formation: one from bud mesenchyme to odontoblasts during embryonic development, and another from apical papilla to odontoblasts during postnatal and adult stages. To quantitatively characterize these developmental paths, we will employ Slingshot to infer pseudotime orderings along these trajectories.</p>"}
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
