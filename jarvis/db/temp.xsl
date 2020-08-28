<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">

<html> 
<head>
<style>
.centered {
  margin: auto;
  width: 90%;
 border: 2px solid black;
  padding: 1px;
}
.right {
  position: absolute;
  right: 0px;
  width: 300px;
  border: 3px solid #73AD21;
  padding: 10px;
}
.left {
  position: absolute;
  left: 0px;
  width: 300px;
  border: 3px solid #73AD21;
  padding: 10px;

}

</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.4.0/3Dmol-min.js"></script>
<script src="https://cdn.plot.ly/plotly-latest.min.js"> </script>
 


<!-- Bootstrap CSS START-->
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous" />


</head>


<body>


  
<!-- Basic table start -->
  <table class="centered table  btn-table" >

    <xsl:for-each select="basic_info/main_relax_info">
    <tr>
      <td>ID: <xsl:value-of select="id"/></td>
      <td>Functional: <xsl:value-of select="method"/></td>
      <td>Primitive cell</td>
      <td>Primitive cell</td>
      <td>Conventional cell</td>
      <td>Conventional cell</td>
    </tr>
    <tr>
      <td>Chemical formula: <xsl:value-of select="formula"/></td>
      <td>Formation energy/atom (eV): <xsl:value-of select="formation_energy"/></td>
      <td>a <xsl:value-of select="a_conv"/> &#8491;</td>
      <td>&#945;:<xsl:value-of select="alpha_conv"/> &#176;</td>
      <td>a <xsl:value-of select="a_prim"/> &#8491;</td>
      <td>&#945;:<xsl:value-of select="alpha_prim"/> &#176;</td>
    </tr>
    <tr>
      <td>Space-group: <xsl:value-of select="spacegroup_number"/></td>
      <td>Relaxed energy/atom (eV): <xsl:value-of select="relaxed_energy"/></td>
      <td>b <xsl:value-of select="b_conv"/> &#8491;</td>
      <td>&#946;:<xsl:value-of select="beta_conv"/> &#176;</td>
      <td>b <xsl:value-of select="a_prim"/> &#8491;</td>
      <td>&#946;:<xsl:value-of select="alpha_prim"/> &#176;</td>
    </tr>
    <tr>
      <td>Material  type: <xsl:value-of select="material_type"/></td>
      <td>SCF bandgap (eV): <xsl:value-of select="scf_indir_gap"/></td>
      <td>c <xsl:value-of select="c_conv"/> &#8491;</td>
      <td>&#947;:<xsl:value-of select="gamma_conv"/> &#176;</td>
      <td>c <xsl:value-of select="c_prim"/> &#8491;</td>
      <td>&#947;:<xsl:value-of select="gamma_prim"/> &#176;</td>
    </tr>
    <tr>
      <td>Crystal system: <xsl:value-of select="crys_system"/></td>
      <td>Point group: <xsl:value-of select="point_group_symbol"/></td>
      <td>Density (gcm<sup>-3</sup>): <xsl:value-of select="density"/></td>
      <td>Volume (<span>&#8491;</span><sup>3</sup>): <xsl:value-of select="volume"/></td>
      <td>nAtoms_prim:<xsl:value-of select="prim_natoms"/></td>
      <td>nAtoms_conv:<xsl:value-of select="conv_natoms"/></td>
    </tr>    
    
    
    
    
    </xsl:for-each>
  </table>
  
  
<!-- Basic table end -->
<br></br>
<!-- Structure viewer start -->
<h3 style="text-align:center">Visualizing atomic structure</h3>
  <div >

  <div 
  class="centered"
  id="geometry"
  style="height: 400px; width: 400px; position: relative;"
  data-backgroundcolor="0xffffff"
  >  
  <span
    style="
      display: inline-block;
      z-index: 1;
      position: absolute;
      bottom: 0px;
      right: 0px;
      font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    "
    >
    
    </span>
  </div>

  
  

<script>
  var data =
    <xsl:value-of select="basic_info/main_relax_info/xyz"/>;
  var viewer_first = $3Dmol.createViewer("geometry", {
    backgroundColor: "white",
  });
  viewer_first.setStyle({stick:{colorscheme:"Jmol"}});
  viewer_first.addModel(data, "xyz");
  viewer_first.setStyle({ sphere: { radius: 0.7 } });
  viewer_first.addStyle({ stick: { radius: 0.2 } });
  viewer_first.zoomTo();
  viewer_first.setClickable({},true,function(atom,viewer,event,container) {
                       viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
                   });
 
  viewer_first.render();
</script>


  </div>
 <!-- Structure viewer end -->
 <br></br>

 
 


<!--Structure info statrts-->

<script>
  var x= <xsl:value-of select="basic_info/main_relax_info/rdf_hist"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Atomic structure analysis");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
    var divElement = document.createElement("div");
    divElement.id = "structure";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script>



<script>
function structure_plotly(){

var data = [data1, data2];
var plot_font = 14;
var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.1, 0.45],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  yaxis1:{tickfont: {size: plot_font,color:'black'},title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  yaxis2: {tickfont: {size: plot_font,color:'black'},anchor: 'x2',title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  xaxis2: {tickfont: {size: plot_font,color:'black'},domain: [0.55, .99],title:{text: 'K-Point',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  
};

Plotly.newPlot('structure', data, layout_convg,{displaylogo: false});

};

   var rdf_hist= <xsl:value-of select="basic_info/main_relax_info/rdf_hist"/>;
 rdf_hist=rdf_hist.split(',').map(Number);
   var rdf_bins= <xsl:value-of select="basic_info/main_relax_info/rdf_bins"/>;
 rdf_bins=rdf_bins.split(',').map(Number);
 var ang_hist1= <xsl:value-of select="basic_info/main_relax_info/ang1_hist"/>;
 ang_hist1=ang_hist1.split(',').map(Number);
   var ang_bins1= <xsl:value-of select="basic_info/main_relax_info/ang1_bins"/>;
 ang_bins1=ang_bins1.split(',').map(Number);
  
  var data1 = {
    x: rdf_bins ,
    y:rdf_hist,
    xaxis: "x1",
    yaxis: "y1",
    type: 'bar',
    marker: {line:{width:.5}}
  };

   var data2 = {
    x: ang_bins1 ,
    y: ang_hist1,
    xaxis: "x2",
    yaxis: "y2",
    type: 'bar',
   marker: {line:{width:.5}}
  };



  structure_plotly();
</script>

<!--Structure end-->  
<br></br>
  
<!--DOS start-->  

<script>
  var x= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/edos_energies"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Electronic density of states");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
    var divElement = document.createElement("div");
    divElement.id = "dos";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
    
    
</script>







<script>
function dos_plotly(){
//  var data=[];
//var data = [data1,data1a];
//var data = [data1,data1a, data2,data2a, data3,data3a,data3b,data3c];
var plot_font = 14;
var layout_convg = {
  grid: {rows: 1, columns: 3},
  xaxis1: {visible:'legendonly',domain: [0, 0.33],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
 
  xaxis2: {visible:'legendonly',domain: [0.33, 0.66],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  xaxis3: {domain: [0.66, 0.99],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
 
  //title:str1.join(cnvg_kp.toString()).join(str2).join(cnvg_enc.toString()),
  title: 'BSDOS',
  showlegend: false,
  
};

Plotly.newPlot('dos', data, layout_convg);

};
   var data=[]

  var x1= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/edos_energies"/>;
  x1=x1.split(',').map(Number);
  var y1= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/total_edos_up"/>;
  y1=y1.split(',').map(Number);
  var y1a= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/total_edos_down"/>;
  y1a=y1a.split(',').map(Number);
  
  var data1 = {
    x: x1 ,
    y: y1,
    xaxis: "x1",
    yaxis: "y1",
    type: "scatter",
  };

    var data1a = {
    x: x1 ,
    y: y1a,
    xaxis: "x1",
    yaxis: "y1a",
    type: "scatter",
  };
   data.push(data1);
   data.push(data1a);
  
  var x2= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/edos_energies"/>;
  x2=x2.split(',').map(Number);
  var s_up= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_up_s"/>;
  s_up=s_up.split(',').map(Number);
  var s_dn= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_down_s"/>;
  s_dn=s_dn.split(',').map(Number);
  var p_up= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_up_p"/>;
  p_up=p_up.split(',').map(Number);
  var p_dn= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_down_p"/>;
  p_dn=p_dn.split(',').map(Number);
 
  var d_up= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_up_d"/>;
  d_up=d_up.split(',').map(Number);
  var d_dn= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/spdf_dos/spin_down_d"/>;
  d_dn=d_dn.split(',').map(Number);
 var keys=["spin_up_s","spin_down_s","spin_up_p","spin_down_p","spin_up_d","spin_down_d"];


  var data_tmp = {
    x: x2,
    y: s_up,
    xaxis: "x2",
    yaxis: "y2",
    name:"s_up",
    marker:{color:"blue"},
    type: "scatter",
  };
  data.push(data_tmp);
  var data_tmp = {
    x: x2,
    y: s_dn,
     
    xaxis: "x2",
    yaxis: "y2",
    marker:{color:"blue"},
    name:"s_down",
    type: "scatter",
  };
  data.push(data_tmp);
  var data_tmp = {
    x: x2,
    y: p_up,
    xaxis: "x2",
    yaxis: "y2",
    type: "scatter",
    marker:{color:"red"},
    name:"p_up",
  };
  data.push(data_tmp);
  var data_tmp = {
    x: x2,
    y: p_dn,
    xaxis: "x2",
    yaxis: "y2",
    type: "scatter",
    marker:{color:"red"},
    name:"p_down",
  };
    data.push(data_tmp);
    var data_tmp = {
    x: x2,
    y: d_up,
    xaxis: "x2",
    yaxis: "y2",
    type: "scatter",
   marker:{color:"green"},
    name:"d_up",
  };
  data.push(data_tmp);
  var data_tmp = {
    x: x2,
    y: d_dn,
    xaxis: "x2",
    yaxis: "y2",
    type: "scatter",
       marker:{color:"green"},
    name:"d_down",
  };
  
  data.push(data_tmp);

  var x3= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/edos_energies"/>;
  x3=x3.split(',').map(Number);
  var y3= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/elemental_dos/spin_up_info"/>;
  y3=y3.split(';');
  var y3a= <xsl:value-of select="basic_info/main_relax_info/main_relax_dos/elemental_dos/spin_down_info"/>;
  y3a=y3a.split(';');

 for (var i=0;i&lt;y3.length-1;i++) {
  
  var data3 = {
    x: x3,
    y: y3[i].split('_')[1].split(',').map(Number),
    xaxis: "x3",
    yaxis: "y3",
    type: "scatter",
    name: y3[i].split('_')[0]
  };


 data.push(data3);
};

 for (var i=0;i&lt;y3.length-1;i++) {
  var data3 = {
    x: x3,
    y: y3a[i].split('_')[1].split(',').map(Number),
    xaxis: "x3",
    yaxis: "y3",
    type: "scatter",
    name: y3a[i].split('_')[0]
  };


 data.push(data3);

};

  
  dos_plotly();
</script>
  
<!--DOS end-->  
  <br></br>




<!--BANDS start--> 
<script>
  var x= <xsl:value-of select="basic_info/main_band/main_bands_info/spin_up_bands_x"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Electronic bandstructure");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
        var divElement = document.createElement("div");
    divElement.id = "bands";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    
    };
</script>



<script src="https://cdn.plot.ly/plotly-latest.min.js"> </script>



<script>
 function bandplot(data){
var plot_font = 14;

var layout_convg = {

  xaxis1: { tickmode: "array",tickvals:kp_labels_points,ticktext:kp_labels ,
  domain: [0.2, 0.7],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},

    yaxis1: {
     range: [-4, 4],
     autorange: false,
   },

  showlegend: false,
  
};



//Plotly.newPlot('bands', data,layout_convg);
};

  var x= <xsl:value-of select="basic_info/main_band/main_bands_info/spin_up_bands_x"/>;

  
  
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_band/main_bands_info/spin_up_bands_y"/>;
  y=y.split(';');


  var data=[]
 
  for (var i=0;i&lt;y.length;i++) {
    if (i==0){
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
    legendgroup:"Spin-up",
    name:"Spin-up",
    marker:{color:'blue'},
    
    type: 'scatter',
  };
    }else{
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
    marker:{color:'blue'},
    showlegend: false,
    
   
  
    type: 'scatter',
  };};
    data.push(data1)
  };
  
   var x= <xsl:value-of select="basic_info/main_band/main_bands_info/spin_down_bands_x"/>;
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_band/main_bands_info/spin_down_bands_y"/>;
  y=y.split(';');

  for (var i=0;i&lt;x.length;i++) {
    if (i==0){
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
   
    legendgroup:"Spin-dn",
    name:"Spin-dn",
    type: 'scatter',
   
    marker:{color:'red'},
  };
    }else{
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
  
    showlegend: false,

    type: 'scatter',
      marker:{color:'red'},
  };};
    data.push(data2);
  };
  
  
  var plot_font = 14;
var kp_labels_points=<xsl:value-of select="basic_info/main_band/main_bands_info/kp_labels_points"/>;
kp_labels_points=kp_labels_points.split(',').map(Number);
var kp_labels=<xsl:value-of select="basic_info/main_band/main_bands_info/kp_labels"/>;
kp_labels=kp_labels.split(',');

  var layout_convg = {

xaxis1: {domain: [0.2, 0.7],tickvals:kp_labels_points,ticktext:kp_labels},
  visible:'legendonly',
    yaxis1: {
     range: [-4, 4],
     autorange: false,
   },
  
  
};

Plotly.newPlot('bands', data,layout_convg);

//bandplot(data);
</script>
  
<!--BANDS end-->  
<br></br>




<!--Spillage start--> 


<script>
  var x= <xsl:value-of select="basic_info/main_spillage_info/soc_bands/spin_up_bands_x"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Spin-orbit coupling spillage");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
    var divElement = document.createElement("div");
    divElement.id = "spillage";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script>



<script src="https://cdn.plot.ly/plotly-latest.min.js"> </script>
<script>
 function bandplot(data){
var plot_font = 14;

var layout_convg = {

  xaxis1: { tickmode: "array",tickvals:kp_labels_points,ticktext:kp_labels ,
  domain: [0.2, 0.7],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},

    yaxis1: {
     range: [-4, 4],
     autorange: false,
   },

  showlegend: false,
  
};



//Plotly.newPlot('bands', data,layout_convg);
};

  var x= <xsl:value-of select="basic_info/main_spillage_info/soc_bands/spin_up_bands_x"/>;
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_spillage_info/soc_bands/spin_up_bands_y"/>;
  y=y.split(';');


  var data=[]
 
  for (var i=0;i&lt;y.length;i++) {
    if (i==0){
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
        xaxis: "x1",
    yaxis: "y1",
    legendgroup:"Spin-up",
    name:"Spin-up",
    marker:{color:'blue'},
    
    type: 'scatter',
  };
    }else{
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
    marker:{color:'blue'},
    showlegend: false,
    
       xaxis: "x1",
    yaxis: "y1",
  
    type: 'scatter',
  };};
    data.push(data1)
  };
  
   var x= <xsl:value-of select="basic_info/main_spillage_info/soc_bands/spin_up_bands_x"/>;
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_spillage_info/soc_bands/spin_up_bands_y"/>;
  y=y.split(';');

  for (var i=0;i&lt;x.length;i++) {
    if (i==0){
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
   
    legendgroup:"Spin-dn",
    name:"Spin-dn",
    type: 'scatter',
       xaxis: "x1",
    yaxis: "y1",
    marker:{color:'red'},
  };
    }else{
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
  
    showlegend: false,
    xaxis: "x1",
    yaxis: "y1",
    type: 'scatter',
      marker:{color:'red'},
  };};
    data.push(data2);
  };
  
  
  
  
  
  
  
  
  
  
   var x= <xsl:value-of select="basic_info/main_spillage_info/nonsoc_bands/spin_up_bands_x"/>;
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_spillage_info/nonsoc_bands/spin_up_bands_y"/>;
  y=y.split(';');


  
 
  for (var i=0;i&lt;y.length;i++) {
    if (i==0){
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
        xaxis: "x2",
    yaxis: "y2",
    legendgroup:"Spin-up",
    name:"Spin-up",
    marker:{color:'blue'},
    
    type: 'scatter',
  };
    }else{
    var data1 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
    marker:{color:'blue'},
    showlegend: false,
    
       xaxis: "x2",
    yaxis: "y2",
  
    type: 'scatter',
  };};
    data.push(data1)
  };
  
   var x= <xsl:value-of select="basic_info/main_spillage_info/nonsoc_bands/spin_up_bands_x"/>;
  x=x.split(';');
  
  var y= <xsl:value-of select="basic_info/main_spillage_info/nonsoc_bands/spin_up_bands_y"/>;
  y=y.split(';');

  for (var i=0;i&lt;x.length;i++) {
    if (i==0){
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
   
    legendgroup:"Spin-dn",
    name:"Spin-dn",
    type: 'scatter',
       xaxis: "x2",
    yaxis: "y2",
    marker:{color:'red'},
  };
    }else{
    var data2 = {
    x: x[i].split(',').map(Number),
    y: y[i].split(',').map(Number),
  
    showlegend: false,
    xaxis: "x2",
    yaxis: "y2",
    type: 'scatter',
      marker:{color:'red'},
  };};
    data.push(data2);
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  var plot_font = 14;
var kp_labels_points=<xsl:value-of select="basic_info/main_spillage_info/soc_bands/kp_labels_points"/>;
kp_labels_points=kp_labels_points.split(',').map(Number);
var kp_labels=<xsl:value-of select="basic_info/main_spillage_info/soc_bands/kp_labels"/>;
kp_labels=kp_labels.split(',');

  var layout_convg = {

xaxis1: {domain: [0.0, 0.5],tickvals:kp_labels_points,ticktext:kp_labels},
  visible:'legendonly',
    yaxis1: {
     range: [-4, 4],
     autorange: false,
   },
  
  xaxis2: {domain: [0.5, 1.0],tickvals:kp_labels_points,ticktext:kp_labels},
  visible:'legendonly',
    yaxis2: {
     range: [-4, 4],
     autorange: false,
   },
};

Plotly.newPlot('spillage', data,layout_convg);
//bandplot(data);
</script>
  
<!--Spillage end-->  
<br></br>


































<!--ggaoptics start--> 
 
<script>
  var x= <xsl:value-of select="basic_info/main_optics_semilocal/main_optics_info/energies"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Optoelectric properties (semi-local)");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
        var divElement = document.createElement("div");
    divElement.id = "ggaoptics";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script>

<script>
function ggaoptics_plotly(){

var data = [data1, data2];
var plot_font = 14;
var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.1, 0.45],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  yaxis1:{tickfont: {size: plot_font,color:'black'},title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  yaxis2: {tickfont: {size: plot_font,color:'black'},anchor: 'x2',title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  xaxis2: {tickfont: {size: plot_font,color:'black'},domain: [0.55, .99],title:{text: 'K-Point',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  
};

Plotly.newPlot('ggaoptics', data, layout_convg);

};

  var energies= <xsl:value-of select="basic_info/main_optics_semilocal/main_optics_info/energies"/>;
  energies=energies.split(',').map(Number);
  var epsilon1= <xsl:value-of select="basic_info/main_optics_semilocal/main_optics_info/real_1"/>;
  epsilon1=epsilon1.split(',').map(Number);
  var epsilon2= <xsl:value-of select="basic_info/main_optics_semilocal/main_optics_info/imag_1"/>;
  epsilon2=epsilon2.split(',').map(Number);
  
  var data1 = {
    x: energies ,
    y:epsilon1,
    xaxis: "x1",
    yaxis: "y1",
    type: 'scatter',
    
  };

   var data2 = {
    x: energies ,
    y: epsilon2,
    xaxis: "x2",
    yaxis: "y2",
    type: 'scatter',
   marker: {line:{width:.5}}
  };



  ggaoptics_plotly();
</script>

<!--ggaoptics end-->  
<br></br>



<!--mbjoptics start-->  



<script>
  var x= <xsl:value-of select="basic_info/main_optics_mbj/main_optics_mbj_info/energies"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Optoelectric properties (TBmBJ)");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
            var divElement = document.createElement("div");
    divElement.id = "mbjoptics";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script>




<script>
function mbjoptics_plotly(){

var data = [data1, data2];
var plot_font = 14;
var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.1, 0.45],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  yaxis1:{tickfont: {size: plot_font,color:'black'},title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  yaxis2: {tickfont: {size: plot_font,color:'black'},anchor: 'x2',title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  xaxis2: {tickfont: {size: plot_font,color:'black'},domain: [0.55, .99],title:{text: 'K-Point',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  
};

Plotly.newPlot('mbjoptics', data, layout_convg);

};

  var energies= <xsl:value-of select="basic_info/main_optics_mbj/main_optics_mbj_info/energies"/>;
  energies=energies.split(',').map(Number);
  var epsilon1= <xsl:value-of select="basic_info/main_optics_mbj/main_optics_mbj_info/real_1"/>;
  epsilon1=epsilon1.split(',').map(Number);
  var epsilon2= <xsl:value-of select="basic_info/main_optics_mbj/main_optics_mbj_info/imag_1"/>;
  epsilon2=epsilon2.split(',').map(Number);
  
  var data1 = {
    x: energies ,
    y:epsilon1,
    xaxis: "x1",
    yaxis: "y1",
    type: 'scatter',
    
  };

   var data2 = {
    x: energies ,
    y: epsilon2,
    xaxis: "x2",
    yaxis: "y2",
    type: 'scatter',
   marker: {line:{width:.5}}
  };



  mbjoptics_plotly();
</script>

<!--mbjoptics end--> 
<br></br>


 

<!--DFPT start--> 


<!--DFPT dielectric tensor table starts--> 


<script>
  var x= <xsl:value-of select="basic_info/main_lepsilon_info/epsilon"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("DFPT dielectric tensor");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
         var divElement = document.createElement("table");
    divElement.id = "dfpt_dielectric_tensor";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);   
    
    };
</script>





<script> 
 var table = document.getElementById("dfpt_dielectric_tensor");
 eij=<xsl:value-of select="basic_info/main_lepsilon_info/epsilon"/> 
 eij=eij.split(',');
 var arr=[["j1","j2","j3"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 for (var i=0;i&lt;eij.length;i++) {
 tmp=eij[i].split(';').map(Number);

 arr.push(tmp);
 

 };
     
 var data = [{
  type: 'table',
 header:{values: ["eij","i1","i2","i3"],height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('dfpt_dielectric_tensor',data);
</script> 
<!--DFPT dielectric tensor table ends-->
<br></br>

<!--DFPT piezoelectric tensor table starts--> 



<script>
  var x= <xsl:value-of select="basic_info/main_lepsilon_info/dfpt_piezoelectric_tensor"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("DFPT piezoelectric tensor");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    };
    
    
             var divElement = document.createElement("table");
    divElement.id = "dfpt_piezoelectric_tensor";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement); 
    
    
</script>

<!--<table class="centered" style="border-spacing: 15px" id="dfpt_piezoelectric_tensor">  </table>
-->



<script> 
 var table = document.getElementById("dfpt_piezoelectric_tensor");
 eij=<xsl:value-of select="basic_info/main_lepsilon_info/dfpt_piezoelectric_tensor"/> 
 eij=eij.split(',');
 var arr=[["j1","j2","j3"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 for (var i=0;i&lt;eij.length;i++) {
 tmp=eij[i].split(';').map(Number);

 arr.push(tmp);
 

 };
     
 var data = [{
  type: 'table',
 header:{values: ["eij","i1","i2","i3","i4","i5","i6"],height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('dfpt_piezoelectric_tensor',data);
</script> 
<!--DFPT piezoelectric tensor table ends-->
<br></br>



  
<!--IRPLOT starts-->

  
  <script>
  var x= <xsl:value-of select="basic_info/main_lepsilon_info/ir_intensity"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("DFPT Infrared intensity");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
        var divElement = document.createElement("div");
    divElement.id = "irplot";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script>
  
  


<script>
var ir_data= <xsl:value-of select="basic_info/main_lepsilon_info/ir_intensity"/>;
ir_frequencies=ir_data.split(';')[0];
gga_ir_intensities=ir_data.split(';')[1];

  ir_frequencies=ir_frequencies.split(',').map(Number);
      
  gga_ir_intensities=gga_ir_intensities.split(',').map(Number);
  
  

var data = 
  {
    
    x: ir_frequencies,
    y: gga_ir_intensities,
    type: 'bar',
  };
  var plot_font=14;
  var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.3, 0.75],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  
};

Plotly.newPlot('irplot', [data],layout_convg);
  </script>

<!--IRPLOT ends--> 
<br></br>

<!--Elastic tensor table starts--> 


 
  <script>
  var x= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/cij"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("FD Elastic constant-tensor");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
     var divElement = document.createElement("table");
    divElement.id = "elastic_tensor";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement); 
    
    };
</script> 
<!--<table class="centered" style="border-spacing: 15px" id="elastic_tensor">  </table>
-->



<script> 
 var table = document.getElementById("elastic_tensor");
 elasticij=<xsl:value-of select="basic_info/main_elastic/main_elastic_info/cij"/> 
 elasticij=elasticij.split(';');
 var arr=[["j1","j2","j3","j4","j5","j6"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 for (var i=0;i&lt;elasticij.length;i++) {
 tmp=elasticij[i].split(',').map(Number);

 arr.push(tmp);
 

 };
     
 var data = [{
  type: 'table',
 header:{values: ["Cij","i1","i2","i3","i4","i5","i6"],height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('elastic_tensor',data);
</script> 
<!--Elastic tensor table ends--> 
<br></br>

<!--Elastic tensor phdos and phbstructure starts--> 

 
<script>
  var x= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/phonon_dos_frequencies"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("FD Phonon DOS at gamma-point only");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    var divElement = document.createElement("div");
    divElement.id = "fddos";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);
    
    };
</script> 

 
<script>
var data=[];
function fddos_plotly(data){


var plot_font = 14;
var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.1, 0.45],tickfont: {size: plot_font,color:'black'},title:{text: 'ENCUT (eV)',font:{size: plot_font,color:"black"}}},
  yaxis1:{tickfont: {size: plot_font,color:'black'},title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  yaxis2: {tickfont: {size: plot_font,color:'black'},anchor: 'x2',title:{text: 'Energy (eV)',font:{size: plot_font,color:"black"}}},
  xaxis2: {tickfont: {size: plot_font,color:'black'},domain: [0.55, .99],title:{text: 'K-Point',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  
};

Plotly.newPlot('fddos', data, layout_convg);

};

  var energies= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/phonon_dos_frequencies"/>;
  energies=energies.split(',').map(Number);
  var ph_dos= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/phonon_dos_intensity"/>;
  ph_dos=ph_dos.split(',').map(Number);

  
  var data1 = {
    x: energies ,
    y:ph_dos,
    xaxis: "x1",
    yaxis: "y1",
    type: 'scatter',
    
  };

   var x= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/phonon_bandstructure_distances"/>;
  x=x.split(',').map(Number);
  
  var y2= <xsl:value-of select="basic_info/main_elastic/main_elastic_info/phonon_bandstructure_frequencies"/>;
  y2=y2.split(';');

  
  
var data_tmp = [];
data.push(data1);


    
 
  for (var i=0;i&lt;y2.length;i++) {
    
    var data2 = {
    x: x,
    y: y2[i].split(',').map(Number),
    marker:{color:'blue'},
    showlegend: false,
    mode:'lines',
    xaxis: "x2",
    yaxis: "y2",
  
    type: 'scatter',
  };
    data.push(data2);
  };
  
    

  
  
  
  




  fddos_plotly(data);
</script>
<br></br>

<!--Elastic tensor phdos and phbstructure ends--> 





<!--FD phonons starts 
 <h3 style="text-align:center">FD Phonon modes at gamma-point only</h3> 
<table class="centered" id="fdphon"></table>

<script> 

 eij=<xsl:value-of select="basic_info/fd_phonons"/> ;
 eij=eij.split(',');
 console.log('eij',eij);
 var arr=[];



 for (var i=0;i&lt;eij.length;i++) {


 arr.push(eij[i]);
 

 };
     
 var data = [{
  type: 'table',
 header:{values: ["ranan"],height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: [eij],
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  
  }
}]






var layout = {
   xaxis: {domain: [0.1, 0.5],},
   
  
}

Plotly.newPlot('fdphon',data, layout);
</script> 


FD phonons ends--> 
<!--FD phonon list -->


<!--EFG tensor table starts--> 

 <script>
  var x= <xsl:value-of select="basic_info/efg_raw_tensor"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("EFG tensor");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
             var divElement = document.createElement("table");
    divElement.id = "efg_tensor";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);   
    
    };
</script> 




<script> 
 var table = document.getElementById("efg_tensor");

     
 elasticij=<xsl:value-of select="basic_info/efg_raw_tensor"/>; 
 console.log('efg',elasticij);
 if (elasticij!==''){



 elasticij=elasticij.split(';');
 var mat=[["Element","Wyckoff","xx","xy","xz","yx","yy","yz","zx","zy","zz"],];

 var count=0;
 var tmp=[];
 for (var i=0;i&lt;elasticij.length-1;i++) {
 tmp=elasticij[i].split(',');

 mat.push(tmp);
 

 };

 var data = [{
  type: 'table',
 header:{values: ["EFG"],height: 30,font: {family: "Arial", size: 24, color: [""]}},

  cells: {
    values: (mat),
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 40,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
};

Plotly.newPlot('efg_tensor',data);
};
</script> 
<!--EFG tensor table ends--> 
<br></br>




<!--BolzTrap tensor table starts--> 

 <script>
  var x= <xsl:value-of select="basic_info/main_boltz/boltztrap_info/pseeb"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="centered";
   var text = document.createTextNode("Boltztrap data");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
             var divElement = document.createElement("table");
    divElement.id = "boltz";
    divElement.setAttribute('class', 'centered');
    document.body.appendChild(divElement);   
    
    };
</script> 




<script> 
 var table = document.getElementById("efg_tensor");

     
var pseeb=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/pseeb"/>; 
pseeb=pseeb.split(',').map(Number);
var pcond=seeb=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/pcond"/>;
pcond=pcond.split(',').map(Number);
var ppf=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/ppf"/>;
ppf=ppf.split(',').map(Number);
var pkappa=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/pkappa"/>;
pkappa=pkappa.split(',').map(Number);
 
var nseeb=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/nseeb"/>;
nseeb=nseeb.split(',').map(Number); 
var ncond=seeb=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/ncond"/>;
ncond=ncond.split(',').map(Number);
var npf=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/npf"/>;
npf=npf.split(',').map(Number);
var nkappa=<xsl:value-of select="basic_info/main_boltz/boltztrap_info/nkappa"/>;
nkappa=nkappa.split(',').map(Number);




 elasticij=elasticij.split(';');
 var mat=[["zz","yy","xx"],];



 mat.push(pseeb);
  mat.push(ppf);
   mat.push(pcond);
   
  mat.push(nseeb);
  mat.push(npf);
   mat.push(ncond);
   



 var data = [{
  type: 'table',
 header:{values:["Data","p-Seebeck","p-Power-factor","p-Conductivity","n-Seebeck","n-Power-factor","n-Conductivity"] ,height: 30,font: {family: "Arial", size: 24, color: [""]}},

  cells: {
    values: (mat),
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 40,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
};

Plotly.newPlot('boltz',data);

</script> 
<!--BolzTrap tensor table ends--> 
<br></br>















</body>
</html>
</xsl:template>
</xsl:stylesheet>


