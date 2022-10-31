/* Credits:
Client-side file upload functionality adapted from:
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

		if (e.target.className == "filedrag" || e.target.className == "filedrag_hover") {
			e.stopPropagation();
			e.preventDefault();
			browse_button = e.target.querySelector(".upload_action_button");
			cancel_button = e.target[1];
			if (e.type == "dragover") {
				e.target.className = "filedrag_hover";
				// disable pointer events for browse and cancel buttons
				browse_button.style.pointerEvents = "none";
				cancel_button.style.pointerEvents = "none";
			}
			else {
				e.target.className = "filedrag";
				// re-enable pointer events for browse and cancel buttons
				browse_button.style.pointerEvents = "auto";
				cancel_button.style.pointerEvents = "auto";
			}
		}
	}


	// file selection
	function FileSelectHandler(e) {

		// cancel event and hover styling
		FileDragHover(e);

		// fetch FileList object
		let files = e.target.files || e.dataTransfer.files;
		let filename = get_filename_from_id(e.target.id);

		// process all File objects
		for (var i = 0, f; f = files[i]; i++) {
			UploadFile(f, filename);
		}

		// clear file list to allow user to re-upload certain files
		e.target.value = "";
	}


	// upload files
	function UploadFile(file, filename) {
		var xhr = new XMLHttpRequest();

		filename_base = filename.split('.')[0];

		if (xhr.upload && file.size <= $id("MAX_FILE_SIZE").value) {

			// remove the old progress bar from a previous upload of the same
			// type of file (counts.tsv, coldata.tsv, config.yml or filter.txt)
			var old_file_progress = $id("progress_of_" + filename_base);
			if (old_file_progress) {
				old_file_progress.remove();
			}

			// create progress bar and set id for tracking
			var progress_div_id = "progress_of_" + filename_base + "_div";
			var progress_div = $id(progress_div_id);
			var progress = progress_div.appendChild(document.createElement("div"));
			displayed_filename = file.name
			if (displayed_filename.length > 12) {
				displayed_filename = displayed_filename.substring(0,9) + "...";
			}
			progress.appendChild(document.createTextNode("upload " + displayed_filename));
			progress.id = "progress_of_" + filename_base;

			// hide button and text
			upload_div = $id(filename_base + "_upload_div");
			upload_div.style.display = "none";

			// show progress div and cancel button
			progress_div.style.display = "inline-block";
			$id("cancel_" + filename_base + "_button").style.display = "block";

			// initially set progress bar to grey
			progress_div.style.backgroundPosition = "100% 0";

			// progress bar
			xhr.upload.addEventListener('progress', function(e) {
				var pc = parseInt(100 - (e.loaded / e.total * 100));
				progress_div.style.backgroundPosition = pc + "% 0";
			}, false);

			// file received/failed
			xhr.onreadystatechange = function(e) {
				if (xhr.readyState == 4) {
					if (xhr.status == 200) {
						var response = JSON.parse(xhr.responseText);

						// remove prior status messages
						var file_status_classname = "status_of" + filename_base;
						var old_file_status = document.getElementsByClassName(file_status_classname);
						for (status_message of old_file_status) {
							status_message.remove();
						}

						if (response.error_status) {
							progress.className = "failed";
							Output(
								"<p class=\"" + file_status_classname +
								"\"><strong>Error in " + file.name + ": " +
								response.error_status + "</strong></p>"
							);
						} else {
							progress.className = "success";
							progress.innerHTML = displayed_filename;
						}
						update_next_button();
					}
				}
			};

			// set form data, including user-supplied filename
			user_filename_req_query = "?user_filename=" + file.name;

			// start upload
			upload_url = window.location.origin + UPLOAD_ENDPOINT + user_filename_req_query
			xhr.open("POST", upload_url, true);

			xhr.setRequestHeader("X_FILENAME", filename);
			xhr.send(file);
		}
	}


	// initialize
	function Init() {

		// is XHR2 available?
		var xhr = new XMLHttpRequest();

		var fileselects = document.getElementsByClassName("fileselect");
		var filedrags = document.getElementsByClassName("filedrag");

		for (fileselect of fileselects) {
			fileselect.addEventListener("change", FileSelectHandler, false);
		}

		if (xhr.upload) {
			// file drop
			for (filedrag of filedrags) {
				filedrag.addEventListener("dragover", FileDragHover, false);
				filedrag.addEventListener("dragleave", FileDragHover, false);
				filedrag.addEventListener("drop", FileSelectHandler, false);
				filedrag.style.display = "block";
			}
		}
	}


	function get_filename_from_id(id) {
		if (id.includes("coldata")) {
			return "coldata.tsv";
		}
		if (id.includes("counts")) {
			return "counts.tsv";
		}
		if (id.includes("config")) {
			return "config.yml";
		}
		if (id.includes("filter")) {
			return "filter.txt";
		}
	}


	// call initialization file
	if (window.File && window.FileList && window.FileReader) {
		Init();
	}

})();
