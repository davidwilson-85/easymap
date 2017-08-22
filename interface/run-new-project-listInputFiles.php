<?php

// Set some flag variables needed later. These flgs indicate whether there are files with certain
// file extensions
$flag_generic = false;
$flag_fasta = false;
$flag_fastq = false;

// Set some empty arrays that will contain the names of the files
$files_fasta = array();
$files_fastq = array();
$files_otherFiles = array();
$all_files = array();

// Define arrays that contain the file extensions supported for reference sequence ad read files
$extensions_fasta = array('.fa', '.fas', '.fasta');
$extensions_fastq = array('.fq', 'fasq', '.fastq');

// Get in an array all items in 'user_data' directory excluding references
// to self (.) and parent (..) dirs
$user_input_files = array_slice(scandir('../user_data'), 2);

// Determine whether there are any files at all
if (count($user_input_files) == 0) {
	// Create JSON object to inform that there are no valid files at all 
	$response = json_encode('error(nofiles)');
} else {
	$flag_generic = true;
	
	// Determine whether there are files with the expected file extensions
	foreach ($user_input_files as $filename) {
		// Get the position (0-based) of last dot in file name
		$last_dot_pos = strrpos($filename, '.');
		// Get the extension of the file
		$extension = substr($filename, $last_dot_pos);
		
		// Fin in the differnet file arrays with the approppriate file names
		if (in_array($extension, $extensions_fasta)) {
			$flag_fasta = true;
			array_push($files_fasta, $filename);
		} else if (in_array($extension, $extensions_fastq)) {
			$flag_fastq = true;
			array_push($files_fastq, $filename);
		} else {
			array_push($files_otherFiles, $filename);
		}
	}

	// Create array of arrays and encode it as a JSON object 
	array_push($all_files, $files_fasta, $files_fastq, $files_otherFiles);
	$response = json_encode($all_files);
}

// Send back JSON object
echo $response;