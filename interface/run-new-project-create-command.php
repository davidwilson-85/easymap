<?php

/*
This file creates the command with info from javascript (user interface) and tells the shell to run
easymap.sh (the master .sh program to run easymap workflows)

Warning: If easymap.sh is run from php, the user is www-data (the apache server). Therefore www-data
must have permission to read and write in the appropriate folders.
*/

//$arguments = $_GET['args'];

// Run workflow
//$command = './easymap-tests.sh project-name snp sim se genome.fa n/p n/p n/p n/p chr1.gff n/p 150+e "0,14;1,31;2,33;3,15;4,5;5,2"+1,10000000+r+50 25+200,40+0,0+1+100 n/p n/p n/p oc ref mut f2wt se';
$command = './easymap -n pn -w ins -sim -r mini_chr1 -i pbinprok2.fa -g complete.gff -sm 3 -ss 5+100,0+500,100+1+50+pe';

shell_exec('cd ..; '. $command);



?>
