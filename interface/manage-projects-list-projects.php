<?php

$projects = array_slice(scandir('../user_projects'), 2);

foreach ($projects as $project) {
	$status_file = '../user_projects/'. $project .'/2_logs/status';	
	
	echo '<div style="background-color: rgb(175,247,124); border: solid green 1px; border-radius: 4px; padding: 10px;">';
	
	$status_contents = fopen($status_file, 'r');
	while(!feof($status_contents)) {
		$line = fgets($status_contents);
		
		$line_fields = explode(':', $line);
		if ($line_fields[0] == 'status') {
			$current_status = trim($line_fields[1]);
		}
		
	}
	fclose($status_contents);
	
	// Project name and status
	echo '<h4>Project: '. $project .'</h4>';
	echo 'Project status: '. $current_status .'<br>';
	
	// Folder size
	if ($current_status != 'running') {
		$folder_size_output = shell_exec('du -sh ../user_projects/'. $project);
		$folder_size_output_array = explode('	', $folder_size_output);
		$folder_size = $folder_size_output_array[0];
		echo 'Project size: '. $folder_size .'<br>';
	}
	
	echo '<div class="managing-buttons">';
		
	// Log link
	echo '<a href="view-log.php?p='. $project .'" class="button">View log file</a>';
	
	if ($current_status == 'finished') {
		echo '<a href="view-report.php?p='. $project .'" class="button">View report</a>';
	}
	
	if ($current_status == 'running') {
		echo '<form><input type="button" class="button" onclick="stopProject(\''. $project .'\')" value="Stop execution" /></form>';
	}
	
	if ($current_status != 'running') {
		echo '<form><input type="button" class="button" onclick="removeProject(\''. $project .'\')" value="Remove from disk" /></form>';
	}
	
	echo '</div></div><br>';
}



	
	
//echo $project .' --- Status --- View log file --- View report --- Stop execution --- Remove from disk<br>';
//                     running    +                 -               +                  + > -                
//                     finished   +                 +               -                  + > -
//                     killed     +                 -               -                  + > -
//                     error      +                 -               -                  + > -





?>
