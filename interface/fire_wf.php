<?php

// Define the names of the input files provided by the user
$project_name = 'project';
$workflow = 'ins';
$data_source = 'sim';
$lib_type = 'se';
$ref_seqs = 'gnm.fa';
$ins_seq = 'ins.fa';
$read_s = 'n/p';
$read_f = 'n/p';
$read_r = 'n/p';
$gff_file = 'gff.gff';
$ann_file = 'n/p';
$sim_mut = '1+ins';
$sim_recsel = 'n/p';
$sim_seq = '10+60,25+0,0+1+50';

// Define a string that contais all the arguments
$arguments = $project_name .' '. $workflow .' '. $data_source .' '. $lib_type .' '. 
$ref_seqs .' '. $ins_seq .' '. $read_s .' '. $read_f .' '. $read_r .' '. $gff_file .' '. 
$ann_file .' '. $sim_mut .' '. $sim_recsel .' '. $sim_seq;

// Run workflow
$command = './master.sh '. $arguments;

//echo $command .'<br>';

$output = shell_exec($command);

// Read workflow log file
$loc = $project_name .'/log.log';
$log_contents = fopen($loc, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);

echo 'Done!';

?>
