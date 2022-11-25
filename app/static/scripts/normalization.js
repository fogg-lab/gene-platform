function displayResults() {
    document.getElementById("download-norm-btn").disabled = false;
    document.getElementById("messages").innerHTML = "Normalization complete.";
}

function updateNextButton() {
    // hide submit button if progress counts or coldata not uploaded
    document.getElementById("submit_button").disabled = !(isUploaded("counts") && isUploaded("coldata"));
}

function submit() {
    document.getElementById("download-norm-btn").disabled = true;
    document.getElementById("messages").innerHTML = "Performing normalization...";
    let submitnormReq = new XMLHttpRequest();
    let submit_norm_query = "method=" + document.getElementById("normalization_method").value;
    submitnormReq.open("POST", "/submit-normalization");
    submitnormReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitnormReq.send(submit_norm_query);
    startTaskUpdateRepeater();
    showRuntime();
}
