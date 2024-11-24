selector_to_html = {"a[href=\"#diffusion-map\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Diffusion Map<a class=\"headerlink\" href=\"#diffusion-map\" title=\"Permalink to this heading\">#</a></h2><p>See the notebook for the detailed code : <a class=\"reference internal\" href=\"20241119_traj_reduction_2.html\"><span class=\"doc std std-doc\">Trajectory reduction with Diffusion Map</span></a></p>", "a[href=\"#force-directed-graph-trimap-phate\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Force-directed graph, Trimap, Phate<a class=\"headerlink\" href=\"#force-directed-graph-trimap-phate\" title=\"Permalink to this heading\">#</a></h2><p>See the notebook for the detailed code : <a class=\"reference internal\" href=\"20241118_traj_reduction.html\"><span class=\"doc std std-doc\">Trajectory reduction with Scanpy External API</span></a></p>", "a[href=\"#trajectory-reduction\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory Reduction<a class=\"headerlink\" href=\"#trajectory-reduction\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>During development, cells exist in a continuum of states. Using single-cell RNA sequencing data, we can identify these states and reconstruct the differentiation process through trajectory analysis. However, trajectory inference involves several challenges. For example, how do we determine the start and end points? How can we ensure that pseudotime, which measures differentiation progress, accurately reflects the biological process? In this section, we will employ various methods for trajectory visualization to gain insights into cell differentiation during odontogenesis. This unbiased approach will lay the foundation for subsequent trajectory inference analyses.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>During development, cells exist in a continuum of states. Using single-cell RNA sequencing data, we can identify these states and reconstruct the differentiation process through trajectory analysis. However, trajectory inference involves several challenges. For example, how do we determine the start and end points? How can we ensure that pseudotime, which measures differentiation progress, accurately reflects the biological process? In this section, we will employ various methods for trajectory visualization to gain insights into cell differentiation during odontogenesis. This unbiased approach will lay the foundation for subsequent trajectory inference analyses.</p>"}
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
