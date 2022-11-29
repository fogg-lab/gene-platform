function displayResults() {
    document.getElementById("download-bc-btn").disabled = false;
}

function updateNextButton() {
    // hide next button if progress counts or coldata not uploaded
    document.getElementById("submit_button").disabled = !(isUploaded("counts") && isUploaded("coldata"));
}

function submit() {
    document.getElementById("download-bc-btn").disabled = true;
    let messages = document.getElementById("messages");
    messages.innerHTML = "Performing batch correction...";
    let datatype_select = document.getElementById("data_type")
    let data_type = datatype_select.options[datatype_select.selectedIndex].text
    let reference_level = document.getElementById("reference_level").value;
    let contrast_level = document.getElementById("contrast_level").value;
    let submitbcReq = new XMLHttpRequest();
    let submit_bc_query = `data_type=${data_type}&reference_level=${reference_level}&contrast_level=${contrast_level}&task_id=${getTaskID()}`;
    submitbcReq.open("POST", "/submit-batch-correction");
    submitbcReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitbcReq.send(submit_bc_query);
    startTaskUpdateRepeater();
    showRuntime();
}
