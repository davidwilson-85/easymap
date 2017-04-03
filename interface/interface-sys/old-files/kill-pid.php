<?php


//'pid' stands for process ID in linux and identifies that process.
$pid = $_POST['project'];

// HERE GET PID OF PROJECT

$pid = 8590;

// The kill command terminates the execution of a process. I use the option -TERM to
// kill also all child processes (.sh, bowtie, etc, python scripts...)
$command1 = 'kill -TERM -'. $pid;

echo $command1;

$output1 = shell_exec($command1);


// Besides stopping the process, I have to update the file status of the workflow/project/
// process. The status file contains two attributes: pid, status.
// The attribute status can have 3 values: running, error, killed, finished
//   running: easymap execution has started
//   error: easymap execution has stoped because of an error
//   killed: easymap was stopped by user before completion 
//   finished: easymap execution has completed succesfully

// This bash comand gets the first line of the file 'status' and then rewrites teh file
// with this info and a new line with 'status:killed'. The overall result is that only
// the second line changes from 'status:running' to 'status:killed'
$command2 = 'pid=$(head -n 1 status); printf $pid"\nstatus:killed" > status';

$output2 = shell_exec($command2);

?>
