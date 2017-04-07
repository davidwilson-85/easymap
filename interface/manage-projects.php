<?php

$projects = array_slice(scandir('../user_projects'), 2);

//$what_to_display = array();

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
	
	
	$style_string = 'style="	
							width: 200px;
							padding: 12px 0;
							margin: 0 0;
							border: 1px solid rgb(0,0,255);
							border-radius: 4px;
							/*box-sizing: border-box;*/
							background-color: #3366ff;
							color: white;
							/*float: left;*/
							text-align: center;
						  "';
	
	
	
	
	// Log link
	echo '<div '. $style_string.' ><a href="view-log.php?p='. $project .'">View log file</a></div><br>';
	
	if ($current_status == 'finished') {
		echo '<div '. $style_string.' ><a href="#">View report</a></div>';
	}
	
	if ($current_status == 'running') {
		echo '<div '. $style_string.' onclick="stopProject(\''. $project .'\')">Stop execution</div>'; // Button to kill project
	}
	
	if ($current_status != 'running') {
		echo '<div '. $style_string.' >Remove from disk</div>'; // Button to remove files
	}
	
	echo '</div><br>';
}



	
	
//echo $project .' --- Status --- View log file --- View report --- Stop execution --- Remove from disk<br>';
//                     running    +                 -               +                  + > -                
//                     finished   +                 +               -                  + > -
//                     killed     +                 -               -                  + > -
//                     error      +                 -               -                  + > -





?>
