function updateNextButton() {
    // hide submit button if counts not uploaded
    document.getElementById("submit_button").disabled = !(isUploaded("counts"));
}

function displayResults() {
    // Done computing correlations
    document.getElementById("messages").innerHTML = "";
    getPlots();
}

function getPearsonPlotListener() {
    if (this.responseText) {
        document.getElementById("pearson_plot").innerHTML = `<img src="data:image/png;base64,${this.responseText}"/>`;
    }
}

function getSpearmanPlotListener() {
    if (this.responseText) {
        document.getElementById("spearman_plot").innerHTML = `<img src="data:image/png;base64,${this.responseText}"/>`;
    }
}

function submit() {
    document.getElementById("messages").innerHTML = "Computing sample correlations...";
    corr_select = document.getElementById("corr_method")
    let corr_method = corr_select.options[corr_select.selectedIndex].value;
    let submitCorrReq = new XMLHttpRequest();
    let submit_corr_query = `corr_method=${corr_method}&task_id=${getTaskID()}`;
    submitCorrReq.open("POST", `/submit-rnaseq-correlation`);
    submitCorrReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitCorrReq.send(submit_corr_query);
    startTaskUpdateRepeater();
    showRuntime();
}

function getPlots() {
    // Reset plots
    for (let plot_div of ["pearson_plot", "spearman_plot"]) {
        plot = document.getElementById(plot_div);
        plot.innerHTML = "";
    }

    // Get requested plots
    switch(document.getElementById("corr_method").value) {
        case "pearson":
            var requested_plots = ["pearson"];
            break;
        case "spearman":
            var requested_plots = ["spearman"];
            break;
        case "both":
            var requested_plots = ["pearson", "spearman"];
    }
    for (corr_method of requested_plots) {
        let getPlotsReq = new XMLHttpRequest();
        let getPlotsQuery = `?corr_method=${corr_method}&task_id=${getTaskID()}`;
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
