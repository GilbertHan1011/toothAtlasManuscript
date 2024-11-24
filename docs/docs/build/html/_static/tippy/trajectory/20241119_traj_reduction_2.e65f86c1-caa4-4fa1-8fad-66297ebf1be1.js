selector_to_html = {"a[href=\"#Visualize-components-in-diffusion-map\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Visualize components in diffusion map<a class=\"headerlink\" href=\"#Visualize-components-in-diffusion-map\" title=\"Permalink to this heading\">#</a></h2><p>We created a function to visualized the components in diffusion map.</p>", "a[href=\"#Setting-up-Environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Setting up Environment<a class=\"headerlink\" href=\"#Setting-up-Environment\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Diffusionmap\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Diffusionmap<a class=\"headerlink\" href=\"#Diffusionmap\" title=\"Permalink to this heading\">#</a></h2><p>Diffusion maps were computed using scanpy with two different principal component inputs (10 and 30 PCs). The diffusion components were subsequently embedded into two-dimensional space using force-directed graph visualization.</p>", "a[href=\"#Trajectory-reduction-with-diffusion-map-in-scanpy\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory reduction with diffusion map in scanpy<a class=\"headerlink\" href=\"#Trajectory-reduction-with-diffusion-map-in-scanpy\" title=\"Permalink to this heading\">#</a></h1><h2>Setting up Environment<a class=\"headerlink\" href=\"#Setting-up-Environment\" title=\"Permalink to this heading\">#</a></h2>"}
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
