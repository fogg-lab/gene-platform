// show the rest of the form if data source is selected
update_form();
document.getElementById("data_source").onchange = update_form;


function update_form() {
    let data_source = document.getElementById("data_source").value;
    if (data_source && data_source == "gdc") {
        document.getElementById("gdc_entry").style.display = 'block';
        document.getElementById("gdc_entry_label").style.display = 'block';
        document.getElementById("geo_entry").style.display = 'none';
        document.getElementById("geo_entry_label").style.display = 'none';
    }
    else if (data_source && data_source == "geo") {
        document.getElementById("gdc_entry").style.display = 'none';
        document.getElementById("gdc_entry_label").style.display = 'none';
        document.getElementById("geo_entry").style.display = 'block';
        document.getElementById("geo_entry_label").style.display = 'block';
    }
    if (!data_source) {
        document.getElementById("submit_button").style.display = 'none';
        document.getElementById("submit_button").disabled=true;
        document.getElementById("download-preprocessed-data-btn").style.display = 'none';
        document.getElementById("gdc_entry").style.display = 'none';
        document.getElementById("gdc_entry_label").style.display = 'none';
        document.getElementById("geo_entry").style.display = 'none';
        document.getElementById("geo_entry_label").style.display = 'none';
    }
    else {
        document.getElementById("submit_button").style.display = 'inline-block';
        document.getElementById("submit_button").disabled=false;
        document.getElementById("download-preprocessed-data-btn").style.display = 'inline-block';
    }
}


function preprocessingReqListener() {
    // preprocessing completed
    let status_msg = this.responseText;
    let messages = document.getElementById("messages");
    messages.innerHTML = status_msg;
}


function displayResults() {
    document.getElementById("download-preprocessed-data-btn").disabled = false;
}


function submit() {
    var submitReq = new XMLHttpRequest();
    submitReq.addEventListener("load", submitReqListener);
    submitReq.open("POST", "/confirm-preprocessing-submission");
    submitReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    datasource_select = document.getElementById("data_source");
    let data_source = datasource_select.options[datasource_select.selectedIndex].text;
    let dsets = "";
    if (data_source.toLowerCase().includes("gdc")) {
        dsets = document.getElementById("gdc_entry").value;
    } else if (data_source.toLowerCase().includes("geo")) {
        dsets = document.getElementById("geo_entry").value;
    } else {
        alert("Please select a data source.");
        return;
    }
    let query = `source=${data_source}&dsets=${dsets}&task_id=${getTaskID()}`
    submitReq.send(query);
}


function submissionConfirmed() {
    datasource_select = document.getElementById("data_source");
    let data_source = datasource_select.options[datasource_select.selectedIndex].text;
    let dsets = "";
    if (data_source.toLowerCase().includes("gdc")) {
        dsets = document.getElementById("gdc_entry").value;
    } else if (data_source.toLowerCase().includes("geo")) {
        dsets = document.getElementById("geo_entry").value;
    } else {
        alert("Please select a data source.");
        return;
    }
    let submitPreprocessingReq = new XMLHttpRequest();
    let submit_preprocessing_query = `source=${data_source}&dsets=${dsets}&task_id=${getTaskID()}`;
    submitPreprocessingReq.addEventListener("load", preprocessingReqListener);
    submitPreprocessingReq.open("POST", "/submit-preprocessing");
    submitPreprocessingReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitPreprocessingReq.send(submit_preprocessing_query);
    document.getElementById("download-preprocessed-data-btn").disabled = true;
    startTaskUpdateRepeater();
    showRuntime();
}
