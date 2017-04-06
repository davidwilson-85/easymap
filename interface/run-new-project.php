<?php

/*
This file creates the command with info from javascript (user interface) and tells the shell to run
easymap.sh (the master .sh program to run easymap workflows)

Warning: If easymap.sh is run from fire_wf.php, the user is www-data (the apache server). Therefore www-data
must have permission to read and write in the appropriate folders.
*/

$arguments = $_GET['args'];

// Run workflow
$command = './easymap.sh project-name snp sim se genome.fa n/p n/p n/p n/p chr1.gff n/p 150+e 0,14-1,31-2,33-3,15-4,5-5,2+1,10000000+r+10 2+200,40+0,0+1+100 n/p n/p n/p oc ref mut f2wt se & echo $!';
echo 'OK';

//shell_exec('cd ..; '. $command);



?>
