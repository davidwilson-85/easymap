/*
This .js file simply calls PHP file via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

window.onload = function() {

	var projectName = document.getElementById("projectName").innerHTML;
	
	// Function to communicate html with php via AJAX to read a log file 
	function readReport(projectName) {
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				 document.getElementById("reportInfo").innerHTML = this.responseText;
			}
		};
		xmlhttp.open("GET", "view-report-read.php?p="+projectName, true);
		xmlhttp.send();
	}
	
	readReport(projectName)
}
