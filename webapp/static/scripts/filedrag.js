/* Credits:
filedrag.js - HTML5 File Drag & Drop
Featured on SitePoint.com
Developed by Craig Buckler (@craigbuckler) of OptimalWorks.net
*/

(function() {

	// getElementById
	function $id(id) {
		return document.getElementById(id);
	}


	// output information
	function Output(msg) {
		var m = $id("messages");
		m.innerHTML = msg + m.innerHTML;
	}


	// file drag hover
	function FileDragHover(e) {
		e.stopPropagation();
		e.preventDefault();
		e.target.className = (e.type == "dragover" ? "hover" : "");
	}


	// file selection
	function FileSelectHandler(e) {

		// cancel event and hover styling
		FileDragHover(e);

		// fetch FileList object
		var files = e.target.files || e.dataTransfer.files;

		// process all File objects
		for (var i = 0, f; f = files[i]; i++) {
			UploadFile(f);
		}

		// clear file list to allow user to re-upload certain files
		var file_list = document.getElementById('fileselect');
		if (file_list){
			file_list.value = "";
		}

	}


	// upload files
	function UploadFile(file) {
		var xhr = new XMLHttpRequest();

		xhr.onreadystatechange = function() {
			if (xhr.readyState == 4) {
				if (xhr.status != 200) {
					//error handling code here
				}
				else {
					//TODO: this
				}
			}
		}

		if (xhr.upload && file.size <= $id("MAX_FILE_SIZE").value) {

			// remove the old progress bar from a previous upload of the same
			// type of file (counts.tsv, coldata.tsv, config.yml or filter.txt)
			var standard_filename = get_standardized_filename(file.name);
			var old_file_progress = $id("progress_of_" + standard_filename);
			if (old_file_progress) {
				old_file_progress.remove();
			}

			// create progress bar and set id for tracking
			var progress_div_id = "progress_of_" + standard_filename + "_div"
			var o = $id(progress_div_id);
			var progress = o.appendChild(document.createElement("p"));
			progress.appendChild(document.createTextNode("upload " + file.name));
			progress.id = "progress_of_" + standard_filename;
			
			// show progress div and cancel button
			o.style.display = 'table-cell';
			$id("cancel_" + standard_filename + "_button").style.display = "inline-block";

			// progress bar
			xhr.upload.addEventListener("progress", function(e) {
				var pc = parseInt(100 - (e.loaded / e.total * 100));
				progress.style.backgroundPosition = pc + "% 0";
			}, false);

			// file received/failed
			xhr.onreadystatechange = function(e) {
				if (xhr.readyState == 4) {
					if (xhr.status == 200) {
						var response = JSON.parse(xhr.responseText);
						
						// remove prior status messages
						var file_status_classname = "status_of" + standard_filename
						var old_file_status = document.getElementsByClassName(file_status_classname);
						for (status_message of old_file_status) {
							status_message.remove()
						}

						if (response.error_status) {
							progress.className = "failed";
							o.style.backgroundColor = "#c00";
							Output(
								"<p class=\"" + file_status_classname +
								"\"><strong>Error in " + file.name + ": " +
								response.error_status + "</strong></p>"
							);
						} else {
							progress.className = "success";
							progress.innerHTML = "File uploaded: " + file.name
							o.style.backgroundColor = "#0c0";
						}
						update_next_button();
					}
				}
			};

			// start upload
			xhr.open("POST", $id("upload").action, true);
			xhr.setRequestHeader("X_FILENAME", file.name);
			xhr.send(file);

		}

	}


	// initialize
	function Init() {

		var fileselect = $id("fileselect"),
			filedrag = $id("filedrag"),
			submitbutton = $id("submitbutton");

		// file select
		fileselect.addEventListener("change", FileSelectHandler, false);

		// is XHR2 available?
		var xhr = new XMLHttpRequest();
		if (xhr.upload) {

			// file drop
			filedrag.addEventListener("dragover", FileDragHover, false);
			filedrag.addEventListener("dragleave", FileDragHover, false);
			filedrag.addEventListener("drop", FileSelectHandler, false);
			filedrag.style.display = "block";

			submitbutton.style.display = "none";
		}

	}


	function get_standardized_filename(filename) {
		if (filename.includes("config")) {
			filename = "config";
		} else if (filename.includes("count")) {
			filename = "counts";
		} else if (filename.includes("col")) {
			filename = "coldata";
		} else if (filename.includes('filt')) {
			filename = "filter";
		}
    
    	return filename
	}


	// call initialization file
	if (window.File && window.FileList && window.FileReader) {
		Init();
	}

})();
