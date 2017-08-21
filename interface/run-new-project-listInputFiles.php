<?php

// Set some flag variables needed later. These flgs indicate whether there are files with certain
// file extensions
$flag_generic = false;
$flag_refSeqs = false;
$flag_readFiles = false;

// Define arrays that contain the file extensions supported for reference sequence ad read files
$extensions_refSeqs = array(".fa", ".fas", ".fasta");
$extensions_readFiles = array(".fq", ".fastq");

// Get in an array all items in 'user_data' directory excluding references
// to self (.) and parent (..) dirs
$user_input_files = array_slice(scandir('../user_data'), 2);

// Determine whether there are any files
if f

// Determine whether there are files with the expected file extensions
foreach ($user_input_files as $filename) {
	// Get the position (0-based) of last dot in file name
	$last_dot_pos = strrpos($filename, ".");
	// Get the extension of the file
	$extension = substr($filename, $last_dot_pos);
	// Output <option value> tags with files with the correct file extension
	if (in_array($extension, $extensions_refSeqs)) {
		$flag = true;
	}
}

// If files with the expected extension exist in ../user_data folder, display a menu to select multiples files; if not, 
// display a warning message
if ($flag == true) {
	echo "<select multiple>";
	// For each files with .fa, .fas, or .fasta extension, create an HTML tag to send to the interface
	foreach ($user_input_files as $filename) {
		if (in_array(substr($filename, strrpos($filename, ".")), $extensions_refSeqs)) {
			echo '<option value="'. $filename .'">'. $filename .'</option>';
		}
	}
	echo "</select>";
} else {
	echo '<div class="warningMessage" style="display:block;">There are no files available with the expected extensions (".fa", ".fas", ".fasta"). Please go to page "Manage input files" to fix this.</div>';
}


?>
