selector_to_html = {"a[href=\"#tuning-other-parameters\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Tuning other parameters<a class=\"headerlink\" href=\"#tuning-other-parameters\" title=\"Permalink to this heading\">#</a></h2><p><a class=\"reference internal\" href=\"20241018_hypertune_para.html\"><span class=\"doc std std-doc\">Tuning with scvi-tools</span></a></p>", "a[href=\"#hyperparameter-tuning\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Hyperparameter tuning<a class=\"headerlink\" href=\"#hyperparameter-tuning\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Now that we have known the best integration method, we need to know what parameters we should set for the integration.\nFortunately, <a class=\"reference external\" href=\"https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html\">scvi-tools</a> offers a hyperparameter tuning tool for users to find the best parameters.</p><p>Here, we devided the hyperparameter tuning into two parts:</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Now that we have known the best integration method, we need to know what parameters we should set for the integration.\nFortunately, <a class=\"reference external\" href=\"https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html\">scvi-tools</a> offers a hyperparameter tuning tool for users to find the best parameters.</p><p>Here, we devided the hyperparameter tuning into two parts:</p>", "a[href=\"#tuning-for-hvg-numbers-and-n-latent\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Tuning for HVG numbers and n_latent<a class=\"headerlink\" href=\"#tuning-for-hvg-numbers-and-n-latent\" title=\"Permalink to this heading\">#</a></h2><p>Firstly, we used grid search methods to select the best parameters for HVG numbers and n_latent. The full script is in this <a class=\"reference internal\" href=\"20241018_gridHypertune.html\"><span class=\"doc std std-doc\">page</span></a>.</p><p>Next, we used <a class=\"reference external\" href=\"https://github.com/YosefLab/scib-metrics\">scib-metrics</a> to plot benchamarking results. <a class=\"reference internal\" href=\"20241019_plot_benchmark.html\"><span class=\"doc std std-doc\">Plot benchmarking results</span></a></p>"}
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
