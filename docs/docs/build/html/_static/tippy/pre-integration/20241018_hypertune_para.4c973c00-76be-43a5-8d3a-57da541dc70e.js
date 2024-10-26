selector_to_html = {"a[href=\"#Layer-number\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Layer number<a class=\"headerlink\" href=\"#Layer-number\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Hypertune-2-:-with-scvi-tools\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Hypertune 2 : with scvi-tools<a class=\"headerlink\" href=\"#Hypertune-2-:-with-scvi-tools\" title=\"Permalink to this heading\">#</a></h1><p>In this script, we used scvi-tools to tune the parameters like n_hidden, n_layers, dropout_rate, gene_likelihood and dispersion. Training loss the the metrics that to evaluate the training performance.</p>", "a[href=\"#n_latent\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">n_latent<a class=\"headerlink\" href=\"#n_latent\" title=\"Permalink to this heading\">#</a></h2><p>Though we runned the latent number here, the results evaluated by scib is more reliable.</p>", "a[href=\"#Dropout-rate\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Dropout rate<a class=\"headerlink\" href=\"#Dropout-rate\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Gene-likelihood\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Gene likelihood<a class=\"headerlink\" href=\"#Gene-likelihood\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Iteration\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Iteration<a class=\"headerlink\" href=\"#Iteration\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Dispersion\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Dispersion<a class=\"headerlink\" href=\"#Dispersion\" title=\"Permalink to this heading\">#</a></h2>"}
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
