document.getElementById("data_type").onchange = show_hide_parameters;

function show_hide_parameters() {

  if (!document.getElementById("data_type").value) {
    document.getElementById("form_parameters").style.display = "none"
    document.getElementById("form_parameters").style.visibility = "hidden"
  }

  else {
    document.getElementById("form_parameters").style.display = "block"
    document.getElementById("form_parameters").style.visibility = "visible"
    
    if (document.getElementById("data_type").value == "RNA-Seq") {
      document.getElementById("use_qual_weights_label").style.display = "none"
      document.getElementById("use_qual_weights_label").style.visibility = "hidden"
      document.getElementById("use_qual_weights").setAttribute("type","hidden")
    }
    else {
      document.getElementById("use_qual_weights_label").style.display = "inline"
      document.getElementById("use_qual_weights_label").style.visibility = "visible"
      document.getElementById("use_qual_weights").setAttribute("type","checkbox")
    }
  }
}

function submitReqListener() {
  console.log(this.responseText);
  // analysis is complete - load the results page
  var displayReq = new XMLHttpRequest()
  window.location.href = "/display"
}

function submit() {
  var submitReq = new XMLHttpRequest();
  submitReq.addEventListener("load", submitReqListener);
  submitReq.open("POST", "/submit");
  submitReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
  var parameters_req_query = ""
  var parameters = get_parameters()
  for (param of parameters) {
    parameters_req_query += param[0]
    parameters_req_query += "="
    parameters_req_query += param[1]
    parameters_req_query += "&"
  }
  // remove last "&" character from the query string
  parameters_req_query = parameters_req_query.slice(0, -1)
  submitReq.send(parameters_req_query);
  show_runtime()
}

function get_parameters() {
  var param_labels = ["padj_thresh", "min_prop", "min_expr", "adj_method", "condition", "contrast_level", "reference_level", "use_qual_weights"]
  var params = []
  datatype_select = document.getElementById("data_type")
  var selection = datatype_select.options[datatype_select.selectedIndex].text
  datatype_keyval = ["data_type", selection]
  params.push(datatype_keyval)

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

function show_runtime() {
  counter = 0
  overlay_div = document.getElementById('runtime_overlay')
  overlay_div.style.display = "flex";
  
  var timer = setInterval(function () {
      document.getElementById("time").innerHTML = "Please wait... Analysis runtime: " + counter + " seconds"
      counter += 1  
  }, 1000);  
}
