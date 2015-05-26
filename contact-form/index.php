<?php
	// Start session.
	session_start();
	
	// Set a key, checked in mailer, prevents against spammers trying to hijack the mailer.
	$security_token = $_SESSION['security_token'] = uniqid(rand());
	
	if ( ! isset($_SESSION['formMessage'])) {
		$_SESSION['formMessage'] = 'Fill in the form below to send me an email.';	
	}
	
	if ( ! isset($_SESSION['formFooter'])) {
		$_SESSION['formFooter'] = ' ';
	}
	
	if ( ! isset($_SESSION['form'])) {
		$_SESSION['form'] = array();
	}
	
	function check($field, $type = '', $value = '') {
		$string = "";
		if (isset($_SESSION['form'][$field])) {
			switch($type) {
				case 'checkbox':
					$string = 'checked="checked"';
					break;
				case 'radio':
					if($_SESSION['form'][$field] === $value) {
						$string = 'checked="checked"';
					}
					break;
				case 'select':
					if($_SESSION['form'][$field] === $value) {
						$string = 'selected="selected"';
					}
					break;
				default:
					$string = $_SESSION['form'][$field];
			}
		}
		return $string;
	}
?><!doctype html>
<html lang="en" class="no-js">
<head>
	
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
		<meta name="generator" content="RapidWeaver" />
		
	<meta charset="UTF-8">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<link href='http://fonts.googleapis.com/css?family=Droid+Sans:400,700|Ubuntu:300,400,700|Merriweather:400,700' rel='stylesheet' type='text/css'>
	<link rel="stylesheet" type="text/css" media="all" href="../rw_common/themes/split/consolidated-12.css" />
		 	
	
	
	<script src="../rw_common/themes/split/js/modernizr.js"></script>
	<title>Contact | Split Navigation for RapidWeaver 6</title>
	
	
	
	
	
</head>

<body>
	<header class="header">
		<div id="logo"><a href="http://demo.yuzoolthemes.com/split/"><img src="../rw_common/images/logo@2x" width="206" height="55" alt="Site logo"/></a></div>

		<nav class="main-nav">
				<ul><li><a href="../" rel="">Home</a></li><li><a href="../styled-3/" rel="" class="parent">About</a></li><li><a href="./" rel="" id="current">Contact</a></li></ul>
		</nav> <!-- main-nav -->
	</header>

	<section id="intro">
		<div id="intro-background"></div>
		<div id="intro-tagline">
			<h1>Split Navigation for RapidWeaver 6</h1> 
		</div> 
	</section> 

	<main class="content">
		<section id="subheader">
		
				<nav class="secondary-nav">
					<ul></ul>
				</nav> <!-- secondary-nav -->
		</section>
		
		<div class="container">
		<aside id="aside" role="complementary">
			<div id="sidebar">
					<h3 class="sidebarTitle"></h3>
					
					<div id="archives">
					
					</div><!-- archives -->
			</div><!-- sidebar -->
		</aside><!-- aside -->
			<div id="contentContainer">
					
<div class="message-text"><?php echo $_SESSION['formMessage']; unset($_SESSION['formMessage']); ?></div><br />

<form action="./files/mailer.php" method="post" enctype="multipart/form-data">
	 <div>
		<label>Your Name</label> *<br />
		<input class="form-input-field" type="text" value="<?php echo check('element0'); ?>" name="form[element0]" size="40"/><br /><br />

		<label>Your Email</label> *<br />
		<input class="form-input-field" type="text" value="<?php echo check('element1'); ?>" name="form[element1]" size="40"/><br /><br />

		<label>Subject</label> *<br />
		<input class="form-input-field" type="text" value="<?php echo check('element2'); ?>" name="form[element2]" size="40"/><br /><br />

		<label>Message</label> *<br />
		<textarea class="form-input-field" name="form[element3]" rows="8" cols="38"><?php echo check('element3'); ?></textarea><br /><br />

		<div style="display: none;">
			<label>Spam Protection: Please don't fill this in:</label>
			<textarea name="comment" rows="1" cols="1"></textarea>
		</div>
		<input type="hidden" name="form_token" value="<?php echo $security_token; ?>" />
		<input class="form-input-button" type="reset" name="resetButton" value="Reset" />
		<input class="form-input-button" type="submit" name="submitButton" value="Submit" />
	</div>
</form>

<br />
<div class="form-footer"><?php echo $_SESSION['formFooter']; unset($_SESSION['formFooter']); ?></div><br />

<?php unset($_SESSION['form']); ?>
			</div>
		</div>
	</main> <!-- content -->
	
	<footer id="footer" role="contentinfo" class="clearfix">	
		<div id="footerContent">
			<div id="footerText">&copy; 2014 RapidWeaver + Yuzoolthemes <a href="#" id="rw_email_contact">Contact Me</a><script type="text/javascript">var _rwObsfuscatedHref0 = "mai";var _rwObsfuscatedHref1 = "lto";var _rwObsfuscatedHref2 = ":he";var _rwObsfuscatedHref3 = "llo";var _rwObsfuscatedHref4 = "@yu";var _rwObsfuscatedHref5 = "zoo";var _rwObsfuscatedHref6 = "lth";var _rwObsfuscatedHref7 = "eme";var _rwObsfuscatedHref8 = "s.c";var _rwObsfuscatedHref9 = "om";var _rwObsfuscatedHref = _rwObsfuscatedHref0+_rwObsfuscatedHref1+_rwObsfuscatedHref2+_rwObsfuscatedHref3+_rwObsfuscatedHref4+_rwObsfuscatedHref5+_rwObsfuscatedHref6+_rwObsfuscatedHref7+_rwObsfuscatedHref8+_rwObsfuscatedHref9; document.getElementById("rw_email_contact").href = _rwObsfuscatedHref;</script></div>
		</div><!-- footerContainer -->
		<div class="clearer"></div>
	</footer><!-- footer -->
	
	<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.0/jquery.min.js"></script>
	<script src="../rw_common/themes/split/js/yuzoolthemes.js"></script>
	<script type="text/javascript" src="../rw_common/themes/split/javascript.js"></script>

</body>
</html><!-- END html -->