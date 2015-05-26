<?php
//start the session
session_start();

//////////////////////////////////////////////////////
// Begin variables to be written out by RapidWeaver //
//////////////////////////////////////////////////////

//set the return URL
$return_url = "../index.php";

//set the users email address
$email = "user@domain.com";

//array of fields in form. (In the format "field_name" => "field_label")
$form_fields = array(
"element0" => 'Your Name',
"element1" => 'Your Email',
"element2" => 'Subject',
"element3" => 'Message'
);

$required_fields = array("element0", "element1", "element2", "element3");

$mail_from_name 	= "element0";
$mail_from_email 	= "element1";
$mail_subject 		= "element2";

//uses the email address defined above as the from email.
$send_from_users_email = false;

//sets the PHP setting 'sendmail_from' for use on a windows server.
$windows_server = false;

// Set up the error and success messages.
$message_success = 'Thank you, your message has been sent.';
$message_unset_fields = "Fields marked with * are required.";

////////////////////////////////////////////////////
// End variables to be written out by RapidWeaver //
////////////////////////////////////////////////////

// Check key variable from form against session key.
if ( ! isset($_POST['form_token']) || $_POST['form_token'] !== $_SESSION['security_token']) {
	// Set a fixed error message if the keys don't match.
	redirect($return_url, 'We cannot verify that you are trying to send an email from this form. Please try again.');
}

// SPAM checking. If the "comment" form field has been filled out,
// send back to form asking to remove content and exit the script.
if ($_POST['comment']) {
	redirect($return_url, 'Please remove content from the last textarea before submitting the form again. This is to protect against SPAM abuse.');
}

/////////////////////////
// PROCESS FORM FIELDS //
/////////////////////////

$magic_quotes = (bool) get_magic_quotes_gpc();
foreach ($_POST['form'] as $key => $value) {
	if ($magic_quotes) {
		$value = stripslashes($value);
	}
	$_SESSION['form'][$key] = $value;
}


///////////////////////////
// CHECK REQUIRED FIELDS //
///////////////////////////

//if any of the required fields are empty
if (check_required_fields($required_fields) === false) {
	//return to form with error message.
	redirect($return_url, $message_unset_fields);
} else {
	
	///////////////////////////////////
	// ALL IS OK, SETUP GLOBAL VAR'S //
	///////////////////////////////////
	
	//check email address
	if ( ! check_email($email)) unset($email);
	
	//set mime boundry. Needed to send the email. Mixed seperates text from attachments.
	$mixed_mime_boundary = 'rms-mix-x'.md5(mt_rand()).'x';
	//alt seperates html from plain text.
	$alt_mime_boundary = 'rms-alt-x'.md5(mt_rand()).'x';
	
	//set the from address if user supplied email is invalid use form owners.
	$submitted_email = '';
	if (isset($_SESSION['form'][$mail_from_email])) {
		$submitted_email = $_SESSION['form'][$mail_from_email];
	}
	
	if (check_email($submitted_email) && $send_from_users_email === false) {
		$from = $reply_to = $_SESSION['form'][$mail_from_name].' <'.$submitted_email.'>';
	} else {
		$from = '<'.$email.'>';
		$reply_to = check_email($submitted_email) ? '<'.$submitted_email.'>' : $from;
	}
	
	//set the email subject
	$subject = '';
	if (isset($_SESSION['form'][$mail_subject])) {
		$subject = $_SESSION['form'][$mail_subject];
	}
	
	//email headers
	if ($windows_server === true) {
		$headers  = "From: $from\r\n" .
		"Reply-to: $reply_to\r\n" .
		"MIME-Version: 1.0\r\nContent-Type: multipart/mixed; " .
		"boundary=$mixed_mime_boundary";
	} else {
		$headers  = "From: $from\n" .
		"Reply-to: $reply_to\n" .
		"MIME-Version: 1.0\nContent-Type: multipart/mixed; " .
		"boundary=$mixed_mime_boundary";
	}
	
	
	
	////////////////////////////
	// CONSTRUCT HTML CONTENT //
	////////////////////////////
	
	//Construct HTML email content, looping through each form element
	
	//Note: When you get to a file attachment you need to use $_FILES['form_element']['name']
	//This will just output the name of the file. The files will actually be attached at the end of the message.
	
	//Set a variable for the message content
	$html_content = "<html>\n<head>\n<title>" .
	safe_escape_string($subject) .
	"</title>\n</head>\n<body>\n<p>\n";
	
	////////////////////////////
	// CONSTRUCT TEXT CONTENT //
	////////////////////////////
	
	//construct a plain text version of the email.
	$text_content = '';
	
	//build a message from the reply for both HTML and text in one loop.
	foreach ($form_fields as $field => $label) {
		$html_content .= '<b>' . safe_escape_string($label) . '</b> ';
		$text_content .= "$label ";
		if (isset($_FILES[$field])) {
			$string = (isset($_FILES[$field]['name'])) ? $_FILES[$field]['name'] : '';
		} else {
			$string = (isset($_SESSION['form'][$field])) ? $_SESSION['form'][$field] : '';
		}
		$html_content .= nl2br(safe_escape_string($string)) . "<br /><br />\n";
		$text_content .= "$string\n\n";
	}
	
	//close the HTML content.
	$html_content .= "</p>\n</body>\n</html>";
	
	/////////////////////////////
	// CONSTRUCT EMAIL MESSAGE //
	/////////////////////////////
	
	//Now we combine both HTML and plain text version of the email into one.
	//Creating the message body which contains a Plain text version and an HTML version,
	//users email client will decide which version to display
	$message = "\r\n--$mixed_mime_boundary\r\n" .
	"Content-Type: multipart/alternative; boundary=$alt_mime_boundary\r\n\r\n" .
	"--$alt_mime_boundary\r\n" .
	"Content-Type: text/plain; charset=UTF-8; format=flowed\r\n" .
	"Content-Transfer-Encoding: Quoted-printable\r\n\r\n" .
	"$text_content\r\n\r\n" .
	"--$alt_mime_boundary\r\n" .
	"Content-Type: text/html; charset=UTF-8\r\n" .
	"Content-Transfer-Encoding: Quoted-printable\r\n\r\n" .
	"$html_content\r\n\r\n" .
	"--$alt_mime_boundary--\r\n\r\n" .
	"\r\n\r\n--$mixed_mime_boundary";
	
	
	//////////////////////
	// FILE ATTACHMENTS //
	//////////////////////
	
	//loop through the $_FILES global array and add each attachment to the form.
	if (isset($_FILES)) {
		foreach ($_FILES as $attachment) {
			$filename = $attachment['name'];
			
			//if the file has been uploaded
			if ($attachment['error'] === UPLOAD_ERR_OK && is_uploaded_file($attachment['tmp_name'])) {
				$file = fopen($attachment['tmp_name'],'rb');
				$data = fread($file,filesize($attachment['tmp_name']));
				fclose($file);
				$data = chunk_split(base64_encode($data));
				
				$message .= "\r\nContent-Type: application/octet-stream; name=\"$filename\"" .
				"\r\nContent-Disposition: attachment; filename=\"$filename\"" .
				"\r\nContent-Transfer-Encoding: base64\r\n\r\n$data\r\n\r\n--$mixed_mime_boundary";
			} else if ($attachment['error'] !== UPLOAD_ERR_NO_FILE) {
				//try to provide a useful error message determined from the error code.
				switch ($attachment['error']) {
					case UPLOAD_ERR_INI_SIZE:
					case UPLOAD_ERR_FORM_SIZE:
						$error = "File $filename exceeds the " . ini_get('upload_max_filesize') . 'B limit for the server.';
						break;
					case UPLOAD_ERR_PARTIAL:
						$error = "Only part of the file $filename could be uploaded, please try again.";
						break;
					default:
						$error = "There has been an error attaching the file $filename, please try again.";
				}
				redirect($return_url, $error);
			}
		}
	}
	
	
	//finish off message
	$message .= '--';
	
	//for windows users.
	if ($windows_server === true) {
		ini_set('sendmail_from', $email);
	}
	
	//if the mail sending works
	if (@mail($email, $subject, $message, $headers)) {
		//set the success message
		$notice = $message_success;
		unset($_SESSION['form']);
	} else {
		$notice = "I'm sorry, there seems to have been an error trying to send your email. Please try again.";
	}
	
	//redirect to the form
	redirect($return_url, $notice);
}

//////////////////////
// GLOBAL FUNCTIONS //
//////////////////////

// Redirects back to the form, an optional session message can be defined.
function redirect($url, $message = NULL) {
	if ($message) {
		$_SESSION['formMessage'] = $message;
	}
	header('Location: ' . $url);
	exit;
}

// Function to escape data inputted from users. This is to protect against embedding
// of malicious code being inserted into the HTML email.
// Sample code: safe_escape_string($_POST['form_field'])
function safe_escape_string($string) {
	return htmlspecialchars($string, ENT_QUOTES);
}

// Function to check the validity of email address.
function check_email($email) {
	return (bool) preg_match('/^([a-z0-9_]|\-|\.)+@(([a-z0-9_]|\-)+\.)+[a-z]{2,}$/i', $email);
}

// Function to check the required fields are filled in.
function check_required_fields($required_fields) {
	foreach ($required_fields as $field) {
		if (( ! isset($_SESSION['form'][$field]) || empty($_SESSION['form'][$field])) &&
			( ! isset($_FILES[$field]) || empty($_FILES[$field]['name']))) {
			return false;
		}
	}
	return true;
}