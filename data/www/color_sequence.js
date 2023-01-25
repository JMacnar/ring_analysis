var ss_color_scheme = {}
ss_color_scheme['-'] = "#606060";
ss_color_scheme['_'] = "#606060";
ss_color_scheme['H'] = "#FF0000";
ss_color_scheme['E'] = "#0000FF";
ss_color_scheme['C'] = "#404040";

var aa_color_scheme = {}
aa_color_scheme['-'] = "#606060";
aa_color_scheme['_'] = "#606060";
aa_color_scheme['M'] = "blue";
aa_color_scheme['L'] = "blue";
aa_color_scheme['V'] = "blue";
aa_color_scheme['I'] = "blue";
aa_color_scheme['W'] = "blue";
aa_color_scheme['A'] = "blue";
aa_color_scheme['F'] = "blue";
aa_color_scheme['T'] = "green";
aa_color_scheme['S'] = "green";
aa_color_scheme['N'] = "green";
aa_color_scheme['Q'] = "green";
aa_color_scheme['H'] = "cyan";
aa_color_scheme['Y'] = "cyan";
aa_color_scheme['E'] = "magenta";
aa_color_scheme['D'] = "magenta";
aa_color_scheme['K'] = "red";
aa_color_scheme['R'] = "red";
aa_color_scheme['G'] = "orange";
aa_color_scheme['P'] = "yellow";

function rc(hash,params) {

  var req = new XMLHttpRequest();
  req.open( "GET", "http://localhost:5000/threading1D?key="+hash+"&params="+params, false );
  req.send( null );
  var e = document.getElementById(hash);
  e.innerHTML = req.responseText;  

  var elements = e.getElementsByClassName('ss-sequence');
  for (var i = 0; i < elements.length; ++i) {
      var item = elements[i]
      item.innerHTML = color_text(item.innerHTML, 'ss');
  }
  elements = e.getElementsByClassName('aa-sequence');
  for (var i = 0; i < elements.length; ++i) {
      var item = elements[i]
      item.innerHTML = color_text(item.innerHTML, 'aa');
  }
}

function color_text(text,which_scale) {

  if(which_scale=='aa') the_scale = aa_color_scheme
  if(which_scale=='ss') the_scale = ss_color_scheme
  out = ""
  for(var j=0; j<text.length; ++j) {
    var c = text[j];
    if(c in the_scale) out = out +"<span style='color:"+the_scale[c]+"'>" + c + "</span>";
    else out = out + c;
  }
  return out;
}

function color_sequence() {

  var elements = document.getElementsByClassName('aa-sequence');
  for (var i = 0; i < elements.length; ++i) {
      item = elements[i];
      item.innerHTML = color_text(item.innerHTML,'aa');
  }

  var elements = document.getElementsByClassName('ss-sequence');
  for (var i = 0; i < elements.length; ++i) {
      item = elements[i];
      item.innerHTML = color_text(item.innerHTML,'ss');
  }
}

