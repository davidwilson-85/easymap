/*
This .js file simply calls PHP files via AJAX
manage-data-files.htm --> manage-data-files.js --> manage-data-files.php --> BASH commands 

This file contains the function to update input file info on screen and to delete files

*/

// Enclose all code inside a onload condition?????????????????????????


// Function to communicate html with php via AJAX to retrieve files info 
function filesInfo() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("filesInfo").innerHTML = this.responseText;
		}
	};
	xmlhttp.open("GET", "manage-input-files.php", true);
	xmlhttp.send();
}

// Call filesInfo() when page loads
filesInfo();



// Call projectsInfo() 3 times separated by one second. This is done everytime page is
// loaded. Before implementing this, I observed that sometimes, when user was redirected
// from run-new-project.htm just after clicking on 'Run new project', the new project
// did not appear in the list. With the current setting, this is solved.
/*
var x = 0;
var intervalID = setInterval(function() {
	projectsInfo();
	if (++x === 3) {
		window.clearInterval(intervalID);
	}
}, 1000);
*/

// Call function filesInfo() every 50 seconds to update projects status regularly
//setInterval(filesInfo, 50000);


// Function to communicate html with php via AJAX to remove a project from disk 
function removeFile(fileName) {
	if (confirm('Are you sure you want to permanently remove this file from disk?')) {
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				// Update files info on screen
				filesInfo()
			}
		};
		xmlhttp.open("GET", "remove-file.php?f="+fileName, true);
		xmlhttp.send();
	}
}
