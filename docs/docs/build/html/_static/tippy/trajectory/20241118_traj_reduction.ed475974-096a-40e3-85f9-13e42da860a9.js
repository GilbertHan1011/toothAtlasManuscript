selector_to_html = {"a[href=\"#Force-directed-graph\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Force-directed graph<a class=\"headerlink\" href=\"#Force-directed-graph\" title=\"Permalink to this heading\">#</a></h2><p>We utilized the <code class=\"docutils literal notranslate\"><span class=\"pre\">draw_graph</span></code> function from the <code class=\"docutils literal notranslate\"><span class=\"pre\">scanpy</span></code> package to generate a force-directed graph layout. Compared to tSNE and UMAP, force-directed graphs provide better representation of global topology, making it our preferred method for visualization</p>", "a[href=\"#Trajectory-Reduction-with-Scanpy-External-API\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory Reduction with Scanpy External API<a class=\"headerlink\" href=\"#Trajectory-Reduction-with-Scanpy-External-API\" title=\"Permalink to this heading\">#</a></h1><h2>Set environment<a class=\"headerlink\" href=\"#Set-environment\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Set-environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Set environment<a class=\"headerlink\" href=\"#Set-environment\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Preprocess\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Preprocess<a class=\"headerlink\" href=\"#Preprocess\" title=\"Permalink to this heading\">#</a></h2><p>We excluded low quality cells here to avoid interference/noise in the downstream analysis</p>", "a[href=\"#phate\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">phate<a class=\"headerlink\" href=\"#phate\" title=\"Permalink to this heading\">#</a></h2><p>We applied <a class=\"reference external\" href=\"https://scanpy.readthedocs.io/en/stable/external/generated/scanpy.external.tl.phate.html\">PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)</a>, a dimensionality reduction method specifically designed to visualize biological progressions in single-cell data by embedding high-dimensional data into two or three dimensions.</p>", "a[href=\"#Trimap\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Trimap<a class=\"headerlink\" href=\"#Trimap\" title=\"Permalink to this heading\">#</a></h2><p>We also performed dimensionality reduction using <a class=\"reference external\" href=\"https://scanpy.readthedocs.io/en/stable/external/generated/scanpy.external.tl.trimap.html\">TriMap</a>, a method that leverages triplet constraints to generate low-dimensional embeddings. This analysis was implemented through the external API functions in scanpy.</p>"}
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
