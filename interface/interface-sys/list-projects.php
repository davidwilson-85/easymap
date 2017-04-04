<?php

$projects = array_slice(scandir('../../user_projects'), 2);

//$what_to_display = array();

foreach ($projects as $project) {
	
	$status_file = '../../user_projects/'. $project .'/2_logs/status';
	echo $status_file .'<br>';
	
	
	
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


	echo 'Status: '. $current_status .'<br>';
	echo '<a href="view-log.php?p=../../user_projects/'. $project .'/2_logs/log.log">View log file</a>';
	
	if ($current_status == 'running') {
		echo ' --- Stop execution'; // Button to kill project
	}
	
	if ($current_status != 'running') {
		echo ' --- Remove files'; // Button to remove files
	}
	
	echo '<br><br>';
}



	//echo $project .''. $current_status .'<br>';
	
	//echo $project .' --- Status --- View log file --- View report --- Stop execution --- Remove files<br>';
	//                     running    +                 -               +                  + > -                
	//                     finished   +                 +               -                  + > -
	//                     killed     +                 -               -                  + > -
	//                     error      +                 -               -                  + > -





?>
