$(document).ready(function(){
	$("#pic_source").change(function(){
		showPic(this);
		$("#startpoint, #endpoint").addClass("invisible");
		$("#status").html("<span class='tempoutput'>Status: waiting...</span>");
		$("#submitajax_est").prop("disabled", true);
		$("#coordx1, #coordy1, #coordx2, #coordy2").empty();
	});
	
	$("#picture").load(function() {
		var canvas_fg = $("#canvas_fg");
		var canvas_bg = $("#canvas_bg");
        context_fg = canvas_fg[0].getContext('2d');
		context_bg = canvas_bg[0].getContext('2d');
		
		context_fg.canvas.height = $(this).height();
        context_fg.canvas.width = $(this).width();
		context_bg.canvas.height = $(this).height();
        context_bg.canvas.width = $(this).width();
	});
	
	$("#err_deviation").hide();
	$("#var_dev").change(function(){
		if ($(this).val() == "value")
			$("#err_deviation").show();
		else
			$("#err_deviation").hide();
	});
	
	$("#slicing").change(function(){
		if ($(this).val() == "star") {
			$("#thickness").prop("disabled", false);
		} else {
			$("#thickness").prop("disabled", true);
		}
	});
	
	/*$(".inputparam").change(function(){
		if (!$("#endpoint").hasClass("invisible")) {
			$("#submitajax_est").prop("disabled", true);
		}
	});*/

	$("#canvas_fg").mousedown(function(event) {
		if ($("#picture").attr("src") == "#") {
			return;
		}
		
		var ratio = $("#picture").width() / $("#picture")[0].naturalWidth;
		
		if ($("#coordx1").is(':empty')) {
			$("#coordx1").text(((event.pageX - $(this).offset().left) / ratio).toFixed()) ;
			$("#coordy1").text(((event.pageY - $(this).offset().top) / ratio).toFixed()) ;
			$("#startpoint").removeClass("invisible");
		} else if ($("#coordx2").is(':empty')) {
			$("#coordx2").text(((event.pageX - $(this).offset().left) / ratio).toFixed());
			$("#coordy2").text(((event.pageY - $(this).offset().top) / ratio).toFixed()) ;
			$("#endpoint").removeClass("invisible");
			
			$("#submitajax_est").removeAttr("disabled");
			
			// draw a rectangle
			var x1 = parseInt($("#coordx1").text())*ratio;
			var y1 = parseInt($("#coordy1").text())*ratio;
			var x2 = parseInt($("#coordx2").text())*ratio;
			var y2 = parseInt($("#coordy2").text())*ratio;
			
			//createLine(x1,y1, x2,y2, "LightGreen");
			createLine(x1,y1, x1,y2, "LightGreen", "canvas_bg");
			createLine(x1,y2, x2,y2, "LightGreen", "canvas_bg");
			createLine(x2,y2, x2,y1, "LightGreen", "canvas_bg");
			createLine(x1,y1, x2,y1, "LightGreen", "canvas_bg");
		}
	});
	
	$("#resetbutton").click(function(){
		clearCanvas("canvas_bg");
		clearCanvas("canvas_fg");
		$("#submitajax_est").prop("disabled", true);
		$("#startpoint, #endpoint").addClass("invisible");
		$("#coordx1, #coordy1, #coordx2, #coordy2").empty();
	});
	
	$("#resetoutputbutton").click(function(){
		$("#status").html("<span class='tempoutput'>Status: waiting...</span>");
		//$("#data_plot").empty();
	});
	
	/*$("#submitajax_prepare").click(function(){
		$("#submitajax_est").removeAttr("disabled");
	});*/
		
	$("#submitajax_est").click(function(){
		$(".tempoutput").remove();
		if (!$("#status").is(':empty'))
			$("#status").prepend("<hr>");
		$("#status").prepend("<span class='tempoutput'><img src='pictures/spinner.gif'> Calculation in progress...<br></span>");
		$("#resetoutputbutton").addClass("disabled");
						
		var req = ocpu.call("load.image.area", {
			picture : $("#pic_source")[0].files[0],
			dataBright : $("#data_bright").val(),
			representation : $("#representation").val(),
			error : $("#error").val(),
			var : Math.pow($("#err_deviation").val(), 2),
			varEst : $("#var_dev").val(),
			confLevel : 0.95,
			levelsOfGray : $("#levelsOfGray").val(),
			boxSize : $("#boxSize").val(),
			thickness : $("#thickness").val(),
			nrSlices : $("#nrSlices").val(),
			slicing : $("#slicing").val(),
			parallel : $("#parallel").val(),
			// x and y switched because of R's coordinate system!!!
			startX : $("#coordy1").text(),
			startY : $("#coordx1").text(),
			endX : $("#coordy2").text(),
			endY : $("#coordx2").text()
		}, function(session) {
			session.getObject(function(outtxt){
				$("#status").prepend(outtxt[0]);
				var ratio = $("#picture").width() / $("#picture")[0].naturalWidth;				
				clearCanvas("canvas_fg");				
				
				createEllipse(outtxt[3]*ratio, outtxt[2]*ratio, outtxt[5]*ratio, outtxt[4]*ratio, outtxt[6]*ratio, "Red", "canvas_fg");
				
				var points = outtxt[7].split(', ');		
				for (i = 0; i<points.length; i+=4) {
					x1 = parseInt(points[i+1])*ratio;
					y1 = parseInt(points[ i ])*ratio;
					x2 = parseInt(points[i+3])*ratio;
					y2 = parseInt(points[i+2])*ratio;
					
					//createLine(x1, y1, x2, y2, "Red", "canvas_fg");
					createPoint(x1, y1, 0.5, "DarkRed", "canvas_fg");
					createPoint(x2, y2, 0.5, "DarkRed", "canvas_fg");
				}				
			});
			
			session = session.getKey();
			$(".tempoutput").remove();
			//$("#status").prepend("<br><a target='_blank' href='/ocpu/tmp/"+ session +"/console'>Click for raw R output</a><br>")			
			
			// $("#data_plot").html("<img src='/ocpu/tmp/"+ session +"/graphics/last/png'>");
			
			$("#resetbutton").removeClass("disabled");
			$("#resetoutputbutton").removeClass("disabled");
		});
		
		req.fail(function(){
			clearCanvas("canvas_fg");
			$(".tempoutput").remove();
			$("#status").prepend("R returned an error: " + req.responseText);
			$("#resetbutton").removeClass("disabled");
			$("#resetoutputbutton").removeClass("disabled");
		});
	});
});
