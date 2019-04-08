$(document).ready(function(){
	$("#pic_source").change(function(){
		showPic(this);
		$("#startpoint, #endpoint").addClass("invisible");
		$("#status").html("<span class='tempoutput'>Status: waiting...</span>");
		$("#data_plot").empty();
		// isprazni canvas?
		$("#submitajax_est, #submitajax_test").prop("disabled", true);
		$("#coordx1, #coordy1, #coordx2, #coordy2").empty();
		sessionStorage.removeItem('generatedData');
		sessionStorage.removeItem('dataCenter');
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

	$(".inputparam").change(function(){
		if (!$("#endpoint").hasClass("invisible")) {
			$("#submitajax_est, #submitajax_test").prop("disabled", true);
		}
	});

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

			$("#submitajax_prepare").removeAttr("disabled");

			// draw a line
			var x1 = parseInt($("#coordx1").text())*ratio;
			var y1 = parseInt($("#coordy1").text())*ratio;
			var x2 = parseInt($("#coordx2").text())*ratio;
			var y2 = parseInt($("#coordy2").text())*ratio;

			createLine(x1,y1, x2,y2, "LightGreen", "canvas_bg");
		}
	});

	$("#resetbutton").click(function(){
		clearCanvas("canvas_bg");
		clearCanvas("canvas_fg");
		$("#submitajax_prepare, #submitajax_est, #submitajax_test").prop("disabled", true);
		$("#startpoint, #endpoint").addClass("invisible");
		$("#coordx1, #coordy1, #coordx2, #coordy2").empty();
		sessionStorage.removeItem('generatedData');
		sessionStorage.removeItem('dataCenter');
	});

	$("#resetoutputbutton").click(function(){
		$("#status").html("<span class='tempoutput'>Status: waiting....</span>");
		$("#data_plot").empty();
	});

	$("#submitajax_prepare, #submitajax_est, #submitajax_test").click(function(){
		$(".tempoutput").remove();
		if (!$("#status").is(':empty'))
			$("#status").prepend("<hr>");
		$("#status").prepend("<span class='tempoutput'><img src='pictures/spinner.gif'> Calculation in progress...<br></span>");
		$("#resetoutputbutton").addClass("disabled");

		if (this.id == "submitajax_prepare") {
			testtype = "prepare";
			sessionStorage.removeItem('generatedData');
			sessionStorage.removeItem('dataCenter');
		} else if (this.id == "submitajax_est") {
			testtype = "est";
		} else if (this.id == "submitajax_test") {
			testtype = "test";
		}

		var generateddata = 0, datacenter = 0;
		if (testtype == "est" || testtype == "test") {
			generateddata = sessionStorage.generatedData;
			datacenter = sessionStorage.dataCenter;
		}

		var req = ocpu.call("load.image", {
			generatedData : generateddata,
			dataCenter : datacenter,
			testType : testtype,
			picture : $("#pic_source")[0].files[0],
			dataBright : $("#data_bright").val(),
			error : $("#error").val(),
			var : Math.pow($("#err_deviation").val(), 2),
			varEst : $("#var_dev").val(),
			confLevel : $("#conf_level").val(),
			levelsOfGray : $("#levelsOfGray").val(),
			boxSize : $("#boxSize").val(),
			thickness : $("#thickness").val(),
			// hypothesis testing
			nulla : $("#H0").val(),
			unit : $("#unit").val(),
			alternative : $("#alternative").val(),
			// x and y switched because of R's coordinate system!!!
			startX : $("#coordy1").text(),
			startY : $("#coordx1").text(),
			endX : $("#coordy2").text(),
			endY : $("#coordx2").text()
		}, function(session) {
			session.getObject(function(outtxt){
				$("#status").prepend(outtxt[0]);

				if (testtype == "prepare") {
					$("#submitajax_est, #submitajax_test").removeAttr("disabled");
				}

				var ratio = $("#picture").width() / $("#picture")[0].naturalWidth;

				x1 = parseInt(outtxt[2])*ratio;
				y1 = parseInt(outtxt[1])*ratio;
				x2 = parseInt(outtxt[4])*ratio;
				y2 = parseInt(outtxt[3])*ratio;

				clearCanvas("canvas_fg");
				createLine(x1, y1, x2, y2, "Red", "canvas_fg");

				// Generated data. Must be saved in browser local storage.
				if (!sessionStorage.generatedData) {
					sessionStorage.generatedData = outtxt[5];
					sessionStorage.dataCenter = outtxt[6];
				}
			});

			session = session.getKey();
			$(".tempoutput").remove();
			$("#status").prepend("<br><a target='_blank' href='/ocpu/tmp/"+ session +"/console'>Click for raw R output</a><br>");

			if (testtype == "prepare") {
				// $("#data_plot").html("<img src='/ocpu/tmp/"+ session +"/graphics/last/png'>");
				$("#status").prepend("<img src='/ocpu/tmp/"+ session +"/graphics/last/png'>");
			}

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
