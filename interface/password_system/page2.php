<?php

session_start(); 
$password = $_SESSION['password'];
echo 'Password is: '. $password;

?>

<?php

if (isset($_SESSION['password'])) {

echo 'Session is set.<br>';

} else {

echo 'Session is unset.<br>';

}






if ($password == "elche") { ?>

<html>
<head>
<title>Correct password</title>
</head>
<body>
<p>UMH leaf mutant DB - Correct password</p>
<p><a href="query.php">Go to query.php</a></p>
</body>
</html>

<?php } else { ?>

<html>
<head>
<title>Wrong password</title>
</head>
<body>
<p>UMH leaf mutant DB - Wrong password. <a href="index.htm">Try again</a>.</p>
</body>
</html>

<?php } ?>