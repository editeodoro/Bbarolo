jQuery.noConflict();
jQuery( document ).ready(function($){

	var menu = $('#menu');
	var nav = $('nav.theme');
	var header = $('#headerWrap');
	var menuToggle = 1;
	var notHeader = $('#feature, #container');
	var subUL = nav.find('li').has('ul');
	var subULA = subUL.children('a');

	// Clean Up

	$('h2#sideTitle:empty').remove();

	// Title Spacer
	
	$('#title').css({'padding-right':nav.outerWidth()});

	// Sub page handling

	nav.find('li').has('ul').addClass('parent').on('click','>a',function(){
		
		var thisLink = $(this).attr('href')

		$(document).on("mouseup mouseout",function (e){
			if (!nav.is(e.target)&& nav.has(e.target).length === 0){nav.find('li').removeClass('open closed');}
		});

		// 2nd Click

		if($(this).parent().hasClass('open')){
			window.location = thisLink;
		}else{

		// 1st Click

			// Show Current Tree, if no tree close all li

			if(!$(this).parent('li').parent('ul').parent('li').hasClass('open')){
				nav.find('li').removeClass('open closed');
			}

			// Add open class
			
			$(this).parent().addClass('open');
			nav.find('li').not('.open').addClass('closed');
			return false;
		}
		
	});

	// Nav click funciton
	
	menu.on('click',function(){
		if( menuToggle == 1){
			$('body').addClass('menu');
			header.add(notHeader).add('body').addClass('menu');
			nav.show();
			menuToggle = 2
		}else{
			$('body').removeClass('menu');
			header.add(notHeader).add('body').removeClass('menu');
			nav.hide();
			menuToggle = 1
		}
	});
});