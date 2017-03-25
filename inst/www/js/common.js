function createLine(x1,y1, x2,y2, color, cnv) {
    /*var length = Math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	var angle  = Math.atan2(y2 - y1, x2 - x1) * 180 / Math.PI;
	var transform = 'rotate('+angle+'deg)';

    var line = $('<div>')
        .appendTo('#picdiv')
        .addClass('line')
        .css({
			'position': 'absolute',
			'transform': transform,
			'background': color
        })
        .width(length)
        //.offset({left: x1, top: y1});  // ovo je prije radilo :-/
		.offset({left: Math.min(x1, x2), top: Math.min(y1, y2)});   // fix 29.08.2016.
		
	if (result) {
		line.addClass('resultLine');
	}
		
    return line;*/
	
	var canvas = document.getElementById(cnv);
	var context = canvas.getContext('2d');

	context.beginPath();
	context.moveTo(x1, y1);
	context.lineTo(x2, y2);
	context.lineWidth = 3;
	context.strokeStyle = color;
	context.stroke();
}

function createPoint(x, y, radius, color, cnv) {
	var canvas = document.getElementById(cnv);
	var context = canvas.getContext('2d');
	
	context.beginPath();
	context.arc(x, y, radius, 0, 2*Math.PI, false);
	context.strokeStyle = color;
	context.stroke();
}

function createEllipse(x, y, radiusX, radiusY, rotation, color, cnv) {
	var canvas = document.getElementById(cnv);
	var context = canvas.getContext('2d');

	context.beginPath();
	context.ellipse(x, y, radiusX, radiusY, (-1)*rotation, 0, 2 * Math.PI);
	context.lineWidth = 2;
	context.strokeStyle = color;
	context.stroke();
}

function clearCanvas(cnv) {
	var canvas = document.getElementById(cnv);
	canvas.width = canvas.width;
}

function showPic(input) {
	if (input.files && input.files[0]) {
		var reader = new FileReader();
		reader.onload = function (e) {
			$('#picture').attr('src', e.target.result);
			$("#coordx1, #coordx2").text("");
			$("#submitajax_prepare, #submitajax_est, #submitajax_test").attr("disabled", "disabled");
		}
		reader.readAsDataURL(input.files[0]);
		$("#picture").css("display", "inline");
	}
	
	$("div.line").remove();
	$("#data_plot").text("");
}