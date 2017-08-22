/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH GENERAL FUNCTIONALITY OF THE PAGE, SUSH AS CHECKING CONFIG FILE OR TRIGGERING A PROJECT
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to communicate html with php via AJAX to check if user_projects directory
// is over config/config>user_projects-size-limit
// This function has to be triggered whe page loads
function allowNewProject() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			var response = this.responseText;
			var responseFields = response.split(",");
			if (responseFields[0] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb). Delete past projects to free disk space.";
			} else if (responseFields[0] >= 75 && responseFields[0] < 100) {
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over 75% of the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb).";
			}
			if (responseFields[1] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("simultWarning").style.display = "block";
				document.getElementById("simultWarning").innerHTML = "WARNING: You have reached the maximum number of simultaneously running projects set by the machine administrator (" + responseFields[3] + " projects). Wait until currently running jobs finish.";
			}
		}
	};
	xmlhttp.open("GET", "allow-new-project.php", true);
	xmlhttp.send();
}

// Function to list all genome reference files to be displayed inside <select multiple> tag
// This function has to be triggered whe page loads
function listInputFiles() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {

			// Parse response in JSON format
			var inputFilesresponse = JSON.parse(this.responseText);

			// Create interface to select reference sequence contigs and insertion sequence
			var refSeqsFormField = ['<select multiple id="refSeqsSelector" size=5>'];
			var insSeqFormField = ['<select id="insSeqSelector">'];
			refSeqsFormField.push('<option value="n/p">Select a file</option>');
			insSeqFormField.push('<option value="n/p">Select a file</option>');
			var fasta = inputFilesresponse[0];
			for (i = 0; i < fasta.length; i++) {
				var optionString = '<option value="' + fasta[i] + '">' + fasta[i] + '</option>';
				refSeqsFormField.push(optionString);
				insSeqFormField.push(optionString);
			}
			refSeqsFormField.push('</select>');
			insSeqFormField.push('</select>');
			var refSeqsFormField = refSeqsFormField.join('');
			var insSeqFormField = insSeqFormField.join('');
			document.getElementById("refSeqs").innerHTML = refSeqsFormField;
			document.getElementById("insSeq").innerHTML = insSeqFormField;

			// Create interface to select GFF and ANN files
			var gffFormField = ['<select id="gffSelector">'];
			var annFormField = ['<select id="annSelector">'];
			gffFormField.push('<option value="n/p">Select a file</option>');
			annFormField.push('<option value="n/p">Select a file</option>');
			var otherFiles = inputFilesresponse[2];
			for (i = 0; i < otherFiles.length; i++) {
				var optionString = '<option value="' + otherFiles[i] + '">' + otherFiles[i] + '</option>';
				gffFormField.push(optionString);
				annFormField.push(optionString);
			}
			gffFormField.push('</select>');
			annFormField.push('</select>');
			var gffFormField = gffFormField.join('');
			var annFormField = annFormField.join('');
			document.getElementById("gffFile").innerHTML = gffFormField;
			document.getElementById("annFile").innerHTML = annFormField;
		}
	};
	xmlhttp.open("GET", "run-new-project-listInputFiles.php?args=refSeqs", true);
	xmlhttp.send();
}

// Function to communicate html with php via AJAX to trigger a workflow run 
function runProject() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("commandResponse").innerHTML = this.responseText;
		}
	};
	var arguments = "test-arg";
	xmlhttp.open("GET", "run-new-project-create-command.php?args=" + arguments, true);
	xmlhttp.send();
}

// Trigger required functions when page loads
allowNewProject();
listInputFiles()



///////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH THE DYNAMIC FORMATTING OF THE FORM AND WITH FIELDS VALIDATION
///////////////////////////////////////////////////////////////////////////////////////////////////////


function resetTextField() {
	this.value = "";
}
	
window.onload = function() {
	
	// Functions that need to be declared after document has completely loaded
	
	// Uptdate the command arguments in the screen
	// This function is called in many other functions after updating the value
	// of an arguments.
	function updateCmd() {
			document.getElementById("commandString").innerHTML = cmdArgs;
	}
	
	//
	function verifyProjectName(){
		var text = document.getElementById("form1").projectName.value;
		if(/[^a-zA-Z0-9]/.test( text) ) {
			//alert('Input is not alphanumeric');
			projectNameValidationInfoMessage = 'Input is not alphanumeric';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else if (text == '') {
			//alert('Yoy must give a name to the project');
			projectNameValidationInfoMessage = 'You must give a name to the project';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else {
			cmdArgs[0] = document.getElementById("form1").projectName.value;
			updateCmd()
			projectNameValidationInfoMessage = '';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "none";
		}
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_analysisType() {
		var options = document.getElementsByClassName("analysisType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button1') {
			cmdArgs[1] = 'ins';
			document.getElementById("insSeqField").style.display = "inline";
			document.getElementById("expDataIns").style.display = "inline";
			document.getElementById("expDataSnp").style.display = "none";
			document.getElementById("simDataIns").style.display = "inline";
			document.getElementById("simDataSnp").style.display = "none";
		} else {
			cmdArgs[1] = 'snp';
			document.getElementById("insSeqField").style.display = "none";
			document.getElementById("expDataIns").style.display = "none";
			document.getElementById("expDataSnp").style.display = "inline";
			document.getElementById("simDataIns").style.display = "none";
			document.getElementById("simDataSnp").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_dataSource() {
		var options = document.getElementsByClassName("dataSource");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button3') {
			cmdArgs[2] = 'exp';
			document.getElementById("expDataInterface").style.display = "inline";
			document.getElementById("simDataInterface").style.display = "none";
		} else {
			cmdArgs[2] = 'sim';
			document.getElementById("expDataInterface").style.display = "none";
			document.getElementById("simDataInterface").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_libType() {
		var options = document.getElementsByClassName("libType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button5') {
			cmdArgs[3] = 'se';
			document.getElementById("expDataSingle").style.display = "inline";
			document.getElementById("expDataPaired").style.display = "none";
			document.getElementById("expDataSingleTwosamples").style.display = "inline";
			document.getElementById("expDataPairedTwosamples").style.display = "none";
			document.getElementById("expDataSnpInfo").style.display = "inline";
			document.getElementById("simSeqFL").style.display = "none";
		} else {
			cmdArgs[3] = 'pe';
			document.getElementById("expDataSingle").style.display = "none";
			document.getElementById("expDataPaired").style.display = "inline";
			document.getElementById("expDataSingleTwosamples").style.display = "none";
			document.getElementById("expDataPairedTwosamples").style.display = "inline";
			document.getElementById("expDataSnpInfo").style.display = "inline";
			document.getElementById("simSeqFL").style.display = "block";
		}
		updateCmd();
	}
	
	// Determine all the file names selected, add them to array, and then to argument
	function refSeqs() {
		var contigs = document.getElementById("refSeqsSelector");
		var contigsList = [];
		
		for (var i=0; i<contigs.length; i++) {
			if (contigs[i].selected == true) {
				contigsList.push(contigs[i].value);
			}
		}
		cmdArgs[4] = contigsList;
		updateCmd();
	}
	
	function processSingleSelectors() {
		if (this.id == 'insSeqSelector') {
			cmdArgs[5] = this.value;
		}
		if (this.id == 'gffSelector') {
			cmdArgs[9] = this.value;
		}
		if (this.id == 'annSelector') {
			cmdArgs[10] = this.value;
		}
		if (this.id == 'readsSingle') {
			cmdArgs[6] = this.value;
		}
		if (this.id == 'readsForward') {
			cmdArgs[7] = this.value;
		}
		if (this.id == 'readsReverse') {
			cmdArgs[8] = this.value;
		}
		
		updateCmd();	
	}
	
	// End of functions
	
	
	// Define array with all the command arguments
	var cmdArgs = ['project', 'workflow', 'dataSource', 'libType', 'refSeqs', 'insSeq',
						'readsS', 'readsF', 'readsR', 'gffFile', 'annFile', 'simMut', 'simRecsel',
						'simSeq'];
	
	// React to interactions with text inputs
	// Reset default content when user clicks on input bux
	document.getElementById("form1").projectName.onfocus = resetTextField;
	document.getElementById("form1").simMutNbr.onfocus = resetTextField;
	document.getElementById("form1").simSeqRD.onfocus = resetTextField;
	document.getElementById("form1").simSeqRL.onfocus = resetTextField;
	document.getElementById("form1").simSeqRD.onfocus = resetTextField;
	document.getElementById("form1").simSeqBER.onfocus = resetTextField;
	document.getElementById("form1").simSeqGBS.onfocus = resetTextField;
	
	// Verify input of text fields
	document.getElementById("form1").projectName.onblur = verifyProjectName;
	
	// Create the command string for the first time (for development purposes only)
	updateCmd();
	
	// React to interactions with 2-way selectors
	document.getElementById("button1").onclick = buttons_analysisType;
	document.getElementById("button2").onclick = buttons_analysisType;
	document.getElementById("button3").onclick = buttons_dataSource;
	document.getElementById("button4").onclick = buttons_dataSource;
	document.getElementById("button5").onclick = buttons_libType;
	document.getElementById("button6").onclick = buttons_libType;
	
	// React to interactions with genome contigs selector
	document.getElementById("form1").refSeqsSelector.onblur = refSeqs;
	
	// React to single selectors (insSeq, gffFile, annFile...)
	document.getElementById("form1").insSeqSelector.onblur = processSingleSelectors;
	document.getElementById("form1").gffSelector.onblur = processSingleSelectors;
	document.getElementById("form1").annSelector.onblur = processSingleSelectors;
	
	document.getElementById("form1").readsSingle.onblur = processSingleSelectors;
	document.getElementById("form1").readsForward.onblur = processSingleSelectors;
	document.getElementById("form1").readsReverse.onblur = processSingleSelectors;
	
}

// TO DO:
// Add verification to check that all the three 2-way buttons have been clicked
// When user clicks on submit, check that all validations passed. This could be done with
// a flag variable.
//
//
//
//
//
//
//
//
//
//





