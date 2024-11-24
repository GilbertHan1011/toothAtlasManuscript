selector_to_html = {"a[href=\"#PAGA-Inferecence\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">PAGA Inferecence<a class=\"headerlink\" href=\"#PAGA-Inferecence\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#Connectivity-inference-with-PAGA\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Connectivity inference with PAGA<a class=\"headerlink\" href=\"#Connectivity-inference-with-PAGA\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we employed <a class=\"reference external\" href=\"https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x\">PAGA</a> to analyze cellular connectivity patterns with odontoblasts. This analysis helped identify clusters potentially involved in odontogenesis across different samples.</p>", "a[href=\"#Sample-level-PAGA-inference\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Sample-level PAGA inference<a class=\"headerlink\" href=\"#Sample-level-PAGA-inference\" title=\"Permalink to this heading\">#</a></h2><p>Here, I turned into sample-level PAGA for better understanding differentiation in samples.</p>", "a[href=\"#Environment-Set-Up\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Environment Set Up<a class=\"headerlink\" href=\"#Environment-Set-Up\" title=\"Permalink to this heading\">#</a></h2>"}
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
