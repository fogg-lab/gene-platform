document.getElementById("data_type").onchange = show_hide_parameters;
job_cancelled = false;
counter = 0


function show_hide_parameters() {
  if (!document.getElementById("data_type").value) {
    document.getElementById("form_parameters").style.display = "none";
    document.getElementById("form_parameters").style.visibility = "hidden";
    // disable submit button
    document.getElementById("submit_button").disabled = true;
  }

  else {
    document.getElementById("form_parameters").style.display = "table";
    document.getElementById("form_parameters").style.visibility = "visible";

    if (document.getElementById("data_type").value == "rnaseq") {
      use_qual_weights_label = document.getElementById("use_qual_weights_label");
      if (use_qual_weights_label != null) {
        document.getElementById("use_qual_weights_label").style.display = "none";
        document.getElementById("use_qual_weights_label").style.visibility = "hidden";
        document.getElementById("use_qual_weights").setAttribute("type","hidden");
      }
    }
    else if (document.getElementById("use_qual_weights_label") != null) {
      document.getElementById("use_qual_weights_label").style.display = "inline";
      document.getElementById("use_qual_weights_label").style.visibility = "visible";
      document.getElementById("use_qual_weights").setAttribute("type","checkbox");
    }
    // enable submit button
    document.getElementById("submit_button").disabled = false;
  }
}


function analysisReqListener() {
  status_msg = this.responseText;
  console.log(status_msg);
  // analysis is complete - load the results page
  if (!job_cancelled) {
    window.location.href = "/display";
  }
}


function submissionConfirmed() {
  var analysisReq = new XMLHttpRequest();
  analysisReq.addEventListener("load", analysisReqListener);
  analysisReq.open("POST", "/submit");
  analysisReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");

  data_type = get_data_type();
  data_type_req_query = "data_type=" + data_type;

  analysisReq.send(data_type_req_query);
  showRuntime();
}


function submit() {
  var submitReq = new XMLHttpRequest();
  submitReq.addEventListener("load", submitReqListener);
  submitReq.open("POST", "/confirm-analysis-submission");
  submitReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
  var parameters_req_query = "";
  var parameters = get_parameters();
  for (param of parameters) {
    parameters_req_query += param[0];
    parameters_req_query += "=";
    parameters_req_query += param[1];
    parameters_req_query += "&";
  }
  data_type = get_data_type();
  parameters_req_query += "data_type=" + data_type;
  submitReq.send(parameters_req_query);
}


function get_parameters() {
  var param_labels = ["padj_thresh", "min_prop", "min_expr", "adj_method", "contrast_level", "reference_level", "use_qual_weights"]
  var params = []
  
  for (param_label of param_labels) {
    param_input = document.getElementById(param_label)
    if (param_input != null) {
      if (param_label == "use_qual_weights") {
        var param_keyval = [param_label, param_input.checked]
        params.push(param_keyval)
      }
      else {
        var param_keyval = [param_label, param_input.value]
        params.push(param_keyval)
      }
    }
  }
  
  return params
}


function get_data_type() {
  datatype_select = document.getElementById("data_type")
  var selection = datatype_select.options[datatype_select.selectedIndex].text
  return selection
}
