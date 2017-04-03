<?php


/*

Maybe this is the way:
Make that process-input.sh, simulator.sh, and the worklflows write their PIDs to file
status in the logs folder. Then with PHP I get the PIDs and use something like:

kill -9 <PID>                  (for easymap.sh)
pkill -9 -P <PPID> xxx.sh      (for the workflows) 
pkill -9 -P <PPID> program     (for bowtie, samtools, bcftools, python scripts...)

*/


$pid_easymap = 0;
$pid_simulator = 0;
$pid_workflow = 0;

$project_name = 'pids_project';
$status_file = 'user_projects/'. $project_name .'/2_logs/status';
$status_contents = file_get_contents($status_file);
$status_contents = fopen($status_file, 'r');
while(!feof($status_contents)) {
	$line = fgets($status_contents);
	$line_fields = explode(' ', $line);
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'easymap') {
		$pid_easymap = $line_fields[2];
	}
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'simulator') {
		$pid_simulator = $line_fields[2];
	}
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'workflow') {
		$pid_workflow = $line_fields[2];
	}
}

fclose($status_contents);

echo $pid_easymap .'<br>';
echo $pid_simulator .'<br>';
echo $pid_workflow .'<br>';

$simulator_children = array(
	'sim-mut.py',
	'sim-recsel.py',
	'sim-seq.py'
);

$workflow_children = array(
	'bowtie2-build-s',
	'bowtie2-align-s',
	'samtools',
	'bcftools'
);

foreach ($simulator_children as $simulator_element) {
	echo 'pkill -9 -P '. $pid_simulator .' '. $simulator_element .'<br>';
	shell_exec('pkill -9 -P '. $pid_simulator .' '. $simulator_element);
}

foreach ($workflow_children as $workflow_element) {
	echo 'pkill -9 -P '. $pid_workflow .' '. $workflow_element .'<br>';
	shell_exec('pkill -9 -P '. $pid_workflow .' '. $workflow_element);
}


/*
// example: y -P: match only child processes of the given parent
$command = 'pkill -P <PPID> process-input.sh';


$output1 = shell_exec('pkill -9 easymap.sh; pkill -9 simulator.sh; pkill -9 sim-seq.py; pkill sim-mut.py; pkill sim-recsel.py');


// here change status to killed by appending a line to the bottom: 'status:killed'

*/








?>
