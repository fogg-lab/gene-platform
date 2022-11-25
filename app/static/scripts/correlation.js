function updateNextButton() {
    // hide submit button if counts not uploaded
    document.getElementById("submit_button").disabled = !(isUploaded("counts") && isUploaded("coldata"));
}

function displayResults() {
    // Done computing correlations
    document.getElementById("messages").innerHTML = "";
    getPlots();
}

function getPearsonPlotListener() {
    // Done getting plots
    let plot_div = document.getElementById("pearson_plot");
    plot_div.innerHTML = '<img src="data:image/png;base64,' + this.responseText + '"/>';
}

function getSpearmanPlotListener() {
    // Done getting plots
    let plot_div = document.getElementById("spearman_plot");
    plot_div.innerHTML = '<img src="data:image/png;base64,' + this.responseText + '"/>';
}

function submit() {
    document.getElementById("messages").innerHTML = "Computing sample correlations...";
    corr_select = document.getElementById("corr_method")
    let corr_method = corr_select.options[corr_select.selectedIndex].value;
    let submitCorrReq = new XMLHttpRequest();
    let submit_corr_query = "corr_method=" + corr_method
    submitCorrReq.open("POST", "/submit-rnaseq-correlation");
    submitCorrReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitCorrReq.send(submit_corr_query);
    startTaskUpdateRepeater();
    showRuntime();
}

function getPlots() {
    for (corr_method of ["pearson", "spearman"]) {
        let getPlotsReq = new XMLHttpRequest();
        let getPlotsQuery = "?corr_method=" + corr_method;
        if (corr_method == "pearson") {
            getPlotsReq.addEventListener("load", getPearsonPlotListener);
        } else {
            getPlotsReq.addEventListener("load", getSpearmanPlotListener);
        }
        dest = `/get-correlation-plot${getPlotsQuery}`;
        getPlotsReq.open("GET", dest);
        getPlotsReq.send();
    }
}
