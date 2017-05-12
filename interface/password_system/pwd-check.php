<?php //Code for dealing with password====================================================
session_start();
if (isset($_POST['password'])) {
	$_SESSION['password'] = $_POST['password'];
	$password = $_SESSION['password'];
} else {
	@$password = $_SESSION['password'];
}
?>

<?php if ($password == 'elche') { ?>

<html>
<head>
<title></title>
</head>
<body>
Correct password!<br>
<a href="other-page.php">Other page</a>
</body>
</html>

<?php } else { ?>

<html>
<head>
<title></title>
</head>
<body>
Wrong password. <a href="index.php">Try again</a>
</body>
</html>

<?php } ?>
