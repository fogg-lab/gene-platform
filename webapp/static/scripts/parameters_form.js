document.getElementById("data_type").onchange = show_hide_parameters;

function show_hide_parameters() {

  if (!document.getElementById("data_type").value) {
    document.getElementById("form_parameters").style.display = "none";
    document.getElementById("form_parameters").style.visibility = "hidden";
    // disable submit button
    document.getElementById("submit_button").disabled = true;
  }

  else {
    document.getElementById("form_parameters").style.display = "block";
    document.getElementById("form_parameters").style.visibility = "visible";
    
    if (document.getElementById("data_type").value == "RNA-Seq") {
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
  console.log(this.responseText);
  // analysis is complete - load the results page
  window.location.href = '/display';
}


function submission_confirmed() {
  var analysisReq = new XMLHttpRequest();
  analysisReq.addEventListener("load", analysisReqListener);
  analysisReq.open("POST", "/submit");
  analysisReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");

  data_type = get_data_type();
  data_type_req_query = "data_type=" + data_type;

  analysisReq.send(data_type_req_query);
  show_runtime();
}


function confirm_submission(confirmation_text) {
  overlay = document.getElementById("overlay");
  overlay.style.display="flex";
  dialog_box = document.getElementById("confirm_div_id")
  if (dialog_box == null) {
    dialog_box = document.createElement("div");
  }
  dialog_box.id = "confirm_div_id";
  dialog_box.innerHTML = confirmation_text;
  dialog_box.innerHTML += "<button id='okay_btn_id'>Okay</button>\n";
  dialog_box.innerHTML += "\n<button id='cancel_btn_id'>Cancel</button>\n";
  dialog_box.style.display="block";
  overlay.appendChild(dialog_box)
  document.getElementById('okay_btn_id').onclick = function() {
    dialog_box.style.display="none";
    overlay.style.display="none";
    if (!(confirmation_text.includes("Error"))) {
      submission_confirmed();
    }
  };
  document.getElementById('cancel_btn_id').onclick = function() {
    dialog_box.style.display="none";
    overlay.style.display="none";
  };
}


function submitReqListener() {
  // analysis is complete - load the results page
  var confirmation_text = this.responseText;
  confirm_submission(confirmation_text);
}


function submit() {
  var submitReq = new XMLHttpRequest();
  submitReq.addEventListener("load", submitReqListener);
  submitReq.open("POST", "/confirmsubmission");
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
  var param_labels = ["padj_thresh", "min_prop", "min_expr", "adj_method", "condition", "contrast_level", "reference_level", "use_qual_weights"]
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

function show_runtime() {
  counter = 0
  overlay_div = document.getElementById('overlay')
  overlay_div.style.display = "flex";
  
  var timer = setInterval(function () {
      document.getElementById("time").innerHTML = "Please wait... Analysis runtime: " + counter + " seconds"
      counter += 1
      // 
  }, 1000);  
}
