function updateNextButton() {
    // hide next button if progress counts or coldata not uploaded
    document.getElementById("next_button").disabled = !(isUploaded("counts") && isUploaded("coldata"));
}

function goto(url) {
    window.location.href = url;
}
