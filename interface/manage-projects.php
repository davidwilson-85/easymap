<?php

$projects = array_slice(scandir('../user_projects'), 2);

//$what_to_display = array();

foreach ($projects as $project) {
	
	$status_file = '../user_projects/'. $project .'/2_logs/status';	
	
	$status_contents = fopen($status_file, 'r');
	while(!feof($status_contents)) {
		$line = fgets($status_contents);
		//echo $line .'<br>';
		
		$line_fields = explode(':', $line);
		if ($line_fields[0] == 'status') {
			$current_status = trim($line_fields[1]);
		}
		
	}
	fclose($status_contents);

	echo 'Project: '. $project .'<br>';
	echo 'Status: '. $current_status .' --- ';
	echo '<a href="view-log.php?p='. $project .'">View log file</a><br>';
	
	// Get folder size
	if ($current_status != 'running') {
		$folder_size_output = shell_exec('du -sh ../user_projects/'. $project);
		$folder_size_output_array = explode('	', $folder_size_output);
		$folder_size = $folder_size_output_array[0];
		echo 'Size: '. $folder_size .'<br>';
	}
	
	if ($current_status == 'finished') {
		echo ' --- View report';
	}
	
	if ($current_status == 'running') {
		echo ' --- Stop execution'; // Button to kill project
	}
	
	if ($current_status != 'running') {
		echo ' --- Remove from disk'; // Button to remove files
	}
	
	echo '<br><br>';
}



	//echo $project .''. $current_status .'<br>';
	
	//echo $project .' --- Status --- View log file --- View report --- Stop execution --- Remove from disk<br>';
	//                     running    +                 -               +                  + > -                
	//                     finished   +                 +               -                  + > -
	//                     killed     +                 -               -                  + > -
	//                     error      +                 -               -                  + > -





?>
