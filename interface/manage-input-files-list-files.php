<?php

// Get in an array all items in 'user_data' directory excluding references
// to self (.) and parent (..) dirs
$user_input_files = array_slice(scandir('../user_data'), 2);

// Remove folder 'gnm-ref' from array
// This folder is only useful for the command-line version
foreach ($user_input_files as $key => $item) {
	if(trim($item) == 'gnm_ref') {
		unset($user_input_files[$key]);
	}
}

echo '<hr id="files-separator"></hr>';

foreach ($user_input_files as $user_input_file) {
	
	// Get file size to display it to user
	$file_size_output = shell_exec('du -sh  ../user_data/'. $user_input_file);
	$file_size_output_array = explode('	', $file_size_output);
	$file_size = $file_size_output_array[0];
	
	echo '
		<div id="user-input-files-list">
			<div class="files-item" id="left">
				<h4>'. $user_input_file .'</h4>
			</div>
			<div class="files-item" id="center">
				<h4>'. $file_size .'</h4>
			</div>
			<div class="files-item" id="right">
				<a href="manage-input-files-retrieve-header.php?f='. $user_input_file .'" target="_blank" class="button">Preview</a>
			</div>
			<div class="files-item" id="right">
				<form>
					<input style="margin: 3 0 7px 0" type="button" class="button" onclick="removeFile(\''. $user_input_file .'\')" value="Remove from disk" />
				</form>
			</div>
			<hr id="files-separator"></hr>
		</div>
	';
}




?>
