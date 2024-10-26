selector_to_html = {"a[href=\"#hyperparameter-tuning\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Hyperparameter tuning<a class=\"headerlink\" href=\"#hyperparameter-tuning\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Now that we have known the best integration method, we need to know what parameters we should set for the integration.\nFortunately, <a class=\"reference external\" href=\"https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html\">scvi-tools</a> offers a hyperparameter tuning tool for users to find the best parameters.</p><p>Here, we devided the hyperparameter tuning into two parts:</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Now that we have known the best integration method, we need to know what parameters we should set for the integration.\nFortunately, <a class=\"reference external\" href=\"https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html\">scvi-tools</a> offers a hyperparameter tuning tool for users to find the best parameters.</p><p>Here, we devided the hyperparameter tuning into two parts:</p>", "a[href=\"#id1\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"headerlink\" href=\"#id1\" title=\"Permalink to this heading\">#</a></h2>"}
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
