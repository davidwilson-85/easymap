<?php


/*

Maybe this is the way:
Make that process-input.sh, simulator.sh, and the worklflows write their PIDs to file
status in the logs folder. Then with PHP I get the PIDs and use something like:

kill -9 <PID>                  (for easymap.sh)
pkill -9 -P <PPID> xxx         (for the workflows) 
pkill -9 -P <PPID> program     (for bowtie, samtools, bcftools, python scripts...)

*/




$targets = array(
	'easymap.sh',
	'process-input.sh',
	'fasta-concat.py',
	'verify-input.py',
	'simulator.sh',
	'sim-mut.py',
	'sim-recsel.py',
	'sim-seq.py',
	'workflow-snp', // remove version number from file name
	'workflow-ins', // remove version number from file name
	'bowtie2-build-s',
	'bowtie2-align-s',
	'samtools',
	'bcftools',
);

foreach ($element in $targets) {
	echo $element .'<br>';
}

// example: y -P: match only child processes of the given parent
$command = 'pkill -P <PPID> process-input.sh';


$output1 = shell_exec('pkill easymap.sh; pkill simulator.sh; pkill sim-seq.py; pkill sim-mut.py; pkill sim-recsel.py');


// here change status to killed



?>
