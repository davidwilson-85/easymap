/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

// Function to communicate html with php via AJAX to check if user_projects directory
// is over config/config>user_projects-size-limit
// This function has to be triggered whe page loads
function allowNewProject() {
	//alert('Button pressed');
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			var response = this.responseText;
			//alert(response);
			
			var responseFields = response.split(",");
			
			if (responseFields[0] >= 100) {
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over the maximum space limit set by the machine administrator (" + responseFields[2] + " bytes). Delete past projects to fere disk space.";
			} else if (responseFields[0] >= 75 && responseFields[0] < 100) {
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over 75% of the maximum space limit set by the machine administrator (" + responseFields[2] + " bytes). Wait until currently running jobs finish.";
			}
			
			if (responseFields[1] >= 100) {
				document.getElementById("simultWarning").style.display = "block";
				document.getElementById("simultWarning").innerHTML = "You have reached the maximum number of simultaneously running projects set by the machine administrator (" + responseFields[3] + " projects)";
			}
			
			// Make display:none the div that contains the form (not when >= 75% but < 100)

		}
	};	
	xmlhttp.open("GET", "allow-new-project.php", true);
	xmlhttp.send();
}


// Function to communicate html with php via AJAX to trigger a workflow run 
function runProject() {
	//alert('Button pressed');
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("commandResponse").innerHTML = this.responseText;
		}
	};
	
	var arguments = "test-arg";
	
	xmlhttp.open("GET", "run-new-project.php?args=" + arguments, true);
	xmlhttp.send();
}

// Trigger function when page loads
allowNewProject();
